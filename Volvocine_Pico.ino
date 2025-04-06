#include <Servo.h>
#include <WiFi.h>
#include <WiFiUdp.h>
#include <LittleFS.h>
#include <math.h> 
#include "agent_config.h"
#include "ServerUtils.h"
#include "WiFiManager.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// WiFi設定
const char* ssid = "Buffalo-G-4510";
const char* password = "33354682";

// UDP設定
IPAddress serverIP(192, 168, 13, 98);
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

#define CONTROL_PERIOD_US 2000 // 制御周期 (μs)
#define LOG_BUFFER_SIZE   36000
CompressedLogData logBuffer[LOG_BUFFER_SIZE];
int logIndex = 0;
bool paused = false;
bool lastButtonState = false;

unsigned long prevLoopEndTime = 0;
float phi = 0;
float omega = 3.0f * 3.14f;
float kappa = 1.0f;  // フィードバックゲイン
float alpha = 0.1f;  // 位相遅れ定数
bool bufferOverflowed = false;

// agent_id: 不変なのでRAMで持つだけでOK (送信時にのみ使用)
int agent_id = 0;

// データ保存間隔を設定 (例: 5ループごとに保存)
const int saveInterval = 5;
int loopCounter = 0;

unsigned long lastRequestTime = 0;  // 最後にリクエストを送信した時刻

// ---------------------------------------------------
// 送信バッファをまとめてUDP送信
//   (各パケット先頭に agent_id の1バイトと送信時の時刻4バイトを付加して送る)
// ---------------------------------------------------
void sendLogBuffer() {
  // サーバーが準備できるまで待機
  while (!isServerReady(udp, serverIP, serverPort)) {
    Serial.println("[ERROR] Server not ready. Retrying in 1 second...");
    delay(500);  // 1秒待機
  }

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
  unsigned long dt = now - prevLoopEndTime;
  prevLoopEndTime = now;

  // サーボ制御
  phi += omega * (float)dt / 1e6f;
  float currentSin = cosf(phi);
  myServo.write(110 + 60 * currentSin);

  // データ保存は指定された間隔でのみ実行
  if (loopCounter % saveInterval == 0) {
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
  }

  // 周期制御
  unsigned long elapsed = micros() - now;
  if (elapsed < CONTROL_PERIOD_US) {
    delayMicroseconds(CONTROL_PERIOD_US - elapsed);
  }

  // ループカウンタをインクリメント
  loopCounter++;
}

void setup() {
  pinMode(digitalInputPin, INPUT);
  Serial.begin(115200);
  analogReadResolution(12);
  myServo.attach(22);
  myServo.write(80);

  // WiFi接続
  connectToWiFi(ssid, password);

  udp.begin(12345);

  // サーバー接続先を選択
  IPAddress serverIP1(192, 168, 13, 98);
  IPAddress serverIP2(192, 168, 13, 99);
  
  warmUpUDP(udp, serverIP1, serverPort);  // ServerUtils.cppの関数を呼び出し
  warmUpUDP(udp, serverIP2, serverPort);  
  
  while (true) {
    if (isServerReady(udp, serverIP1, serverPort)) {
      serverIP = serverIP1;
      Serial.println("[INFO] Connected to server at 192.168.13.98");
      break;
    } else if (isServerReady(udp, serverIP2, serverPort)) {
      serverIP = serverIP2;
      Serial.println("[INFO] Connected to server at 192.168.13.99");
      break;
    } else {
      Serial.println("[WARN] No servers are ready. Retrying in 1 second...");
      delay(1000);  // 1秒待機して再試行
    }
  }

  // agent_id 読み込み
  agent_id = readAgentIdFromFile(); // ユーザ実装の想定
  Serial.printf("Loaded agent_id: %d\n", agent_id);

  // サーボモータを真ん中に動かす
  myServo.write(90);
  Serial.println("[INFO] Servo moved to center position (90 degrees)");

  Serial.println("[INFO] Ready to log in RAM");
  prevLoopEndTime = micros();

  // 初期状態をオフに設定
  paused = true;
  logIndex = 0;  // バッファインデックスを初期化
  sendLogBuffer();
  Serial.println("[INFO] System is paused. Press the button to start.");
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

      // サーバーにパラメータをリクエスト
      requestParametersFromServer(udp, serverIP, serverPort, agent_id, omega, kappa, alpha);
      lastRequestTime = millis();  // リクエスト送信時刻を記録
    }
  }
  lastButtonState = currentButtonState;

  // ポーズ中に一定間隔でパラメータをリクエスト
  if (paused && millis() - lastRequestTime >= 5000) {  // 1秒以上経過
    requestParametersFromServer(udp, serverIP, serverPort, agent_id, omega, kappa, alpha);
    lastRequestTime = millis();  // リクエスト送信時刻を更新
  }

  // 記録中
  if (!paused) {
    logSensorData();
  }
}
