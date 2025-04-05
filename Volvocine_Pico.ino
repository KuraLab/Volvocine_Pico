#include <Servo.h>
#include <WiFi.h>
#include <WiFiUdp.h>
#include <LittleFS.h>
#include <math.h>  // for fmodf, M_PI
#include "agent_config.h"

// WiFi設定
const char* ssid = "Buffalo-G-4510";
const char* password = "33354682";

// UDP設定
IPAddress serverIP(192, 168, 13, 5);
unsigned int serverPort = 5000;
WiFiUDP udp;

// ピン設定
const int digitalInputPin = 3;  // ボタン
const int analogPin1 = 27;
const int analogPin2 = 28;

Servo myServo;

// 1レコード5バイトの圧縮構造体 (RAM保持用)
struct __attribute__((packed)) CompressedLogData {
  uint16_t micros16;  // 2バイト: (micros >> 10)
  uint8_t  analog0;   // 1バイト
  uint8_t  analog1;   // 1バイト
  uint8_t  analog2;   // 1バイト
};

#define CONTROL_PERIOD_US 10000
#define LOG_BUFFER_SIZE   36000
CompressedLogData logBuffer[LOG_BUFFER_SIZE];
int logIndex = 0;
bool paused = false;
bool lastButtonState = false;

unsigned long prevLoopEndTime = 0;
float phi = 0;
float omega = 3.0f * 3.14f;
bool bufferOverflowed = false;

// agent_id: 不変なのでRAMで持つだけでOK (送信時にのみ使用)
int agent_id = 0;

// ---------------------------------------------------
// 送信バッファをまとめてUDP送信
//   (各パケット先頭に agent_id の1バイトと送信時の時刻4バイトを付加して送る)
// ---------------------------------------------------
void sendLogBuffer() {
  const int maxPacketBytes = 512;
  uint8_t packet[maxPacketBytes];

  int sentCount = 0;
  int i = 0;

  while (i < logIndex) {
    size_t offset = 0;

    // 1) agent_id (1バイト)
    packet[offset++] = (uint8_t)agent_id;

    // 2) 送信時刻 (4バイト, micros)
    uint32_t sendMicros = micros();
    Serial.printf("[DEBUG] Sending packet at micros=%lu, micros>>10=%lu\n\r", sendMicros, sendMicros >> 10);  // ←追加（送信時刻デバッグ）
    memcpy(&packet[offset], &sendMicros, sizeof(sendMicros));
    offset += sizeof(sendMicros);  // 4バイト

    // 3) ログデータを詰める
    int perPacketCount = 0;
    while (i < logIndex) {
      if (offset + sizeof(CompressedLogData) > maxPacketBytes) {
        break;
      }

      memcpy(&packet[offset], &logBuffer[i], sizeof(CompressedLogData));
      offset += sizeof(CompressedLogData);
      i++;
      perPacketCount++;
    }

    // 4) UDP送信
    udp.beginPacket(serverIP, serverPort);
    udp.write(packet, offset);
    udp.endPacket();

    sentCount += perPacketCount;
  }

  Serial.printf("[INFO] Sent %d records from RAM\n", sentCount);

  if (bufferOverflowed) {
    Serial.println("[WARN] Some data may have been lost due to buffer overflow.");
    bufferOverflowed = false;
  }
}


// ---------------------------------------------------
// センサ読み取り＋RAMバッファ保存（dtはサーボ用のみ）
// ---------------------------------------------------
void logSensorData() {
  unsigned long now = micros();
  Serial.printf("[DEBUG] Sensor logged at micros=%lu, micros>>10=%lu\n\r", now, now >> 10);  // ←追加（センサ取得時刻デバッグ）
  unsigned long dt = now - prevLoopEndTime;
  prevLoopEndTime = now;

  // サーボ制御
  phi += omega * (float)dt / 1e6f;
  float currentSin = cosf(phi);
  myServo.write(110 + 60 * currentSin);

  // ログ用構造体
  CompressedLogData entry;
  entry.micros16 = now >> 10;

  // analog0: phiを [0..2π) → 0..255 に圧縮
  float phiMod = fmodf(phi, 2.0f * (float)M_PI);
  if (phiMod < 0) phiMod += 2.0f * (float)M_PI;
  entry.analog0 = (uint8_t)(phiMod * (255.0f / (2.0f * (float)M_PI)));

  // analog1
  int raw1 = analogRead(analogPin1);  // 0..4095
  entry.analog1 = (uint8_t)(raw1 >> 4);

  // analog2
  int raw2 = analogRead(analogPin2);
  int extended2 = raw2 << 2;  // ×4
  if (extended2 > 4095) extended2 = 4095;
  entry.analog2 = (uint8_t)(extended2 >> 4);

  // バッファに書き込み
  if (logIndex < LOG_BUFFER_SIZE) {
    logBuffer[logIndex++] = entry;
  } else {
    // 1度だけWarnを出す
    if (!bufferOverflowed) {
      Serial.println("[WARN] log buffer overflow!");
      bufferOverflowed = true;
    }
  }

  // バッファ使用率 (10件毎に表示)
  if (logIndex % 10 == 0) {
    float usage = (float)logIndex / LOG_BUFFER_SIZE * 100.0f;
    Serial.printf("[STATUS] buffer: %d/%d (%.1f%%)\n", logIndex, LOG_BUFFER_SIZE, usage);
  }

  // 周期制御
  unsigned long elapsed = micros() - now;
  if (elapsed < CONTROL_PERIOD_US) {
    delayMicroseconds(CONTROL_PERIOD_US - elapsed);
  }
}

void setup() {
  pinMode(digitalInputPin, INPUT);
  Serial.begin(115200);
  analogReadResolution(12);
  myServo.attach(22);

  // WiFi接続
  Serial.print("Connecting to WiFi...");
  WiFi.begin(ssid, password);
  while (WiFi.status() != WL_CONNECTED) {
    delay(500);
    Serial.print(".");
  }
  Serial.println("\nWiFi connected.");
  Serial.print("IP address: ");
  Serial.println(WiFi.localIP());

  // agent_id 読み込み
  agent_id = readAgentIdFromFile(); // ユーザ実装の想定
  Serial.printf("Loaded agent_id: %d\n", agent_id);

  Serial.println("[INFO] Ready to log in RAM");
  prevLoopEndTime = micros();
}

void loop() {
  bool currentButtonState = digitalRead(digitalInputPin);
  if (currentButtonState && !lastButtonState) {
    paused = !paused;
    Serial.println(paused ? "[INFO] Paused - Sending log from RAM" : "[INFO] Resumed");
    delay(300);  // チャタリング防止

    if (paused) {
      // ログ送信
      sendLogBuffer();
      // バッファ初期化
      logIndex = 0;
    }
  }
  lastButtonState = currentButtonState;

  // 記録中
  if (!paused) {
    logSensorData();
  }
}
