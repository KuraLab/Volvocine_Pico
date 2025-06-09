#include <Servo.h>
#include <WiFi.h>
#include <WiFiUdp.h>
#include <LittleFS.h>
#include <math.h> 
#include <vector>
#include <algorithm>
#include <tuple> // std::tupleを使用するために必要
#include <cstdlib>   // rand(), srand()
#include "agent_config.h"
#include "ServerUtils.h"
#include "WiFiManager.h"
#include "calculateTrimmedMean.h"

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

// 1レコード6バイトの圧縮構造体 (RAM保持用)
struct __attribute__((packed)) CompressedLogData {
  uint32_t micros24 : 24;  // 3バイト: (micros >> 8)
  uint8_t  analog0;        // 1バイト
  uint8_t  analog1;        // 1バイト
  uint8_t  analog2;        // 1バイト
};

#define CONTROL_PERIOD_US 2000 // 制御周期 (μs)
#define LOG_BUFFER_SIZE   28000
CompressedLogData logBuffer[LOG_BUFFER_SIZE];
int logIndex = 0;
bool paused = false;
bool lastButtonState = false;

// START受信時のログ開始時刻
unsigned long startLoggingMillis = 0;
unsigned long startLoggingMicros = 0;
float t_delay;

unsigned long prevLoopEndTime = 0;
unsigned long prevLoopEndTime2 = 0;
float phi = 0;
float omega = 3.0f * 3.14f;
float kappa = 1.0f;  // フィードバックゲイン
float kappa_init = 0.0f;
float kappa_now = 0.0f;
float alpha = 0.1f;  // 位相遅れ定数
bool bufferOverflowed = false;
float wait_max = 2.0f * M_PI / omega;


// agent_id: 不変なのでRAMで持つだけでOK (送信時にのみ使用)
int agent_id = 0;

// データ保存間隔を設定 (例: 5ループごとに保存)
const int saveInterval = 5;
int loopCounter = 0;

unsigned long lastRequestTime = 0;  // 最後にリクエストを送信した時刻

// 窓サイズを定義
const int windowSize = 1000; // 必要なサイズに変更
std::vector<int> raw2Window(windowSize, 0); // 固定サイズのリングバッファ
int raw2Index = 0; // 現在のインデックスを管理

// 正規化する関数
float normalize(float value, float lower, float upper) {
  // 正規化
  float normalized = (value - lower) / (upper - lower) - 0.5f;

  // ±0.5にクリップ
  if (normalized > 0.5f) {
    normalized = 0.5f;
  } else if (normalized < -0.5f) {
    normalized = -0.5f;
  }

  return normalized;
}

// ---------------------------------------------------
// 送信バッファをまとめてUDP送信
//   (各パケット先頭に agent_id の1バイトと送信時の時刻4バイトを付加して送る)
// ---------------------------------------------------
void sendLogBuffer() {
  const int maxPacketBytes = 512;
  uint8_t packet[maxPacketBytes];
  const int maxRetries = 100;

  int sentCount = 0;
  int i = 0;

  while (i < logIndex) {
    int retry = 0;
    bool ackReceived = false;

    while (retry < maxRetries && !ackReceived) {
      size_t offset = 0;
      int startIndex = i;

      // サーバー準備チェック
      while (!isServerReady(udp, serverIP, serverPort)) {
        Serial.println("[ERROR] Server not ready. Retrying in 1 second...");
        delay(500);
        if (WiFi.status() != WL_CONNECTED) {
          connectToWiFi(ssid, password);
        }
      }

      // 1) agent_id (1バイト)
      packet[offset++] = (uint8_t)agent_id;

      // 2) 送信時刻 (4バイト)
      uint32_t sendMicros = micros();
      memcpy(&packet[offset], &sendMicros, sizeof(sendMicros));
      offset += sizeof(sendMicros);  // 4バイト

      // 3) データパック詰め
      int perPacketCount = 0;
      uint32_t lastMicros24 = 0;

      while (i < logIndex) {
        if (offset + sizeof(CompressedLogData) > maxPacketBytes) {
          break;
        }

        // タイムスタンプを一時変数経由でコピー
        uint32_t micros24Value = logBuffer[i].micros24;
        memcpy(&packet[offset], &micros24Value, 3);
        offset += 3;

        memcpy(&packet[offset], &logBuffer[i].analog0, sizeof(CompressedLogData) - 3);
        offset += sizeof(CompressedLogData) - 3;

        lastMicros24 = micros24Value;  // 最後の値を保存
        i++;
        perPacketCount++;
      }

      // 4) UDP送信
      udp.beginPacket(serverIP, serverPort);
      udp.write(packet, offset);
      udp.endPacket();

      Serial.printf("[INFO] Packet sent (%d records). Waiting for ACK...\n", perPacketCount);

      // 5) ACK待機
      ackReceived = waitForAck(udp, agent_id, lastMicros24, 1000);
      if (!ackReceived) {
        retry++;
        Serial.printf("[WARN] ACK not received (retry %d/%d). Resending...\n", retry, maxRetries);
        i = startIndex;  // 再送時は戻る
        delay(100);
      } else {
        sentCount += perPacketCount;
      }
    }

    if (!ackReceived) {
      Serial.println("[ERROR] Failed to receive ACK after multiple retries. Aborting this packet.");
    }
  }

  Serial.printf("[INFO] Sent %d records from RAM (with ACK)\n", sentCount);

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
  unsigned long elapsed = now - startLoggingMicros;
  prevLoopEndTime = now;

  int raw1 = analogRead(analogPin1);  // 0..4095
  int raw2 = analogRead(analogPin2);

  // リングバッファにデータを追加
  raw2Window[raw2Index] = raw2;
  raw2Index = (raw2Index + 1) % windowSize; // インデックスを循環させる

  // 下位10%と上位10%の値、およびその平均を計算
  auto [lowerValue, upperValue, trimmedMean] = calculateTrimmedMean(raw2Window, windowSize);

  // 正規化
  float flex = normalize((float)raw2 / 4095.0f, lowerValue, upperValue);

  // サーボ制御
  phi += (kappa_now * cosf((float)elapsed / 1e6f * omega + phi - alpha) * flex) * (float)dt / 1e6f;
  float currentCos = cosf((float)elapsed / 1e6f * omega + phi);
  myServo.write(110 + 60 * currentCos);

  // データ保存は指定された間隔でのみ実行
  if (loopCounter % saveInterval == 0) {
    // ログ用構造体
    CompressedLogData entry;
    entry.micros24 = now >> 8;  // 24ビットに圧縮

    // analog0: phiを [0..2π) → 0..255 に圧縮
    float phiMod = fmodf((float)elapsed / 1e6f * omega + phi, 2.0f * (float)M_PI);
    if (phiMod < 0) phiMod += 2.0f * (float)M_PI;
    entry.analog0 = (uint8_t)(phiMod * (255.0f / (2.0f * (float)M_PI)));

    entry.analog1 = (uint8_t)(raw1 >> 4);

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
      //Serial.printf("[STATUS] buffer: %d/%d (%.1f%%)\n", logIndex, LOG_BUFFER_SIZE, usage);
    }
  }


  unsigned long now2 = micros();
  unsigned long dt2 = now2 - prevLoopEndTime2;
  // 周期制御
  if (dt2 < CONTROL_PERIOD_US) {
    delayMicroseconds(CONTROL_PERIOD_US - dt2);
    prevLoopEndTime2 = micros();
    //Serial.printf("[INFO] Loop took %lu us (expected %d us)\n", dt2, CONTROL_PERIOD_US);
  } else{
    prevLoopEndTime2 = micros();
    //Serial.printf("[WARN] Loop took too long: %lu us (expected %d us)\n", dt2, CONTROL_PERIOD_US);
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

  requestParametersFromServer(udp, serverIP, serverPort, agent_id, omega, kappa, alpha);

  // サーボモータを真ん中に動かす
  myServo.write(90);
  Serial.println("[INFO] Servo moved to center position (90 degrees)");

  Serial.println("[INFO] Ready to log in RAM");
  prevLoopEndTime = micros();
  prevLoopEndTime2 = prevLoopEndTime;

  // 初期状態をオフに設定
  paused = true;
  logIndex = 0;  // バッファインデックスを初期化
  sendLogBuffer();
  kappa_now = kappa_init;
  srand(micros());
  Serial.println("[INFO] System is paused. Press the button to start.");
}

// UDPコマンド受信処理
void checkControlCommand() {
  int packetSize = udp.parsePacket();
  if (packetSize > 0) {
    char buf[16] = {0};
    udp.read(buf, sizeof(buf) - 1);
    if (strcmp(buf, "START") == 0 && paused == true) {
      paused = false;
      startLoggingMillis = millis(); // ログ開始時刻を記録
      startLoggingMicros = micros(); // ログ開始時刻を記録
      Serial.println("[INFO] Received START command from server.");
      t_delay = (rand() / (float)RAND_MAX) * wait_max;
      startLoggingMicros += (unsigned long)(t_delay * 1e6f);
    } else if (strcmp(buf, "STOP") == 0 && paused == false) {
      paused = true;
      Serial.println("[INFO] Received STOP command from server.");
      sendLogBuffer();
      logIndex = 0;
      kappa_now = kappa_init;
    }
  }
}

void loop() {
  checkControlCommand();

  if (!paused && (millis() - startLoggingMillis >= 20000)) {
    kappa_now = kappa;
  } 
  if (!paused && (millis() - startLoggingMillis >= 100000)) {
    while(millis() - startLoggingMillis <= 180000){
      delay(1000);
    }
    paused = true;
    Serial.println("[INFO] Received STOP command from server.");
    sendLogBuffer();
    logIndex = 0;
    kappa_now = kappa_init;
  } 

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
      kappa_now = kappa_init;

      // サーバーにパラメータをリクエスト
      requestParametersFromServer(udp, serverIP, serverPort, agent_id, omega, kappa, alpha);
      lastRequestTime = millis();  // リクエスト送信時刻を記録
    } else{
      startLoggingMillis = millis(); // ログ開始時刻を記録
      startLoggingMicros = micros(); // ログ開始時刻を記録 
      t_delay = (rand() / (float)RAND_MAX) * wait_max;
      startLoggingMicros += (unsigned long)(t_delay * 1e6f);
    }
  }
  lastButtonState = currentButtonState;

  // ポーズ中に一定間隔でパラメータをリクエスト
  if (paused && millis() - lastRequestTime >= 30000) {
    while (WiFi.status() != WL_CONNECTED) {
      connectToWiFi(ssid, password);
    }
    Serial.println("[INFO] WiFi connected.");
    
    requestParametersFromServer(udp, serverIP, serverPort, agent_id, omega, kappa, alpha);
    lastRequestTime = millis();  // リクエスト送信時刻を更新
  }

  // 記録中
  if (!paused) {
    logSensorData();
  }
}
