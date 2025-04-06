#include "ServerUtils.h"
#include <Arduino.h>

bool isServerReady(WiFiUDP& udp, IPAddress serverIP, unsigned int serverPort) {
  const char* handshakeMessage = "HELLO";
  const int timeoutMs = 1000;  // 応答待ちタイムアウト (ミリ秒)
  char response[10];

  // ハンドシェイクメッセージを送信
  udp.beginPacket(serverIP, serverPort);
  udp.write(handshakeMessage);
  udp.endPacket();

  // 応答を待つ
  unsigned long startTime = millis();
  while (millis() - startTime < timeoutMs) {
    int packetSize = udp.parsePacket();
    if (packetSize > 0) {
      int len = udp.read(response, sizeof(response) - 1);
      if (len > 0) {
        response[len] = '\0';  // 文字列終端を追加
        if (strcmp(response, "READY") == 0) {
          Serial.println("[INFO] Server is ready.");
          return true;
        }
      }
    }
  }

  Serial.println("[WARN] No response from server.");
  return false;
}

void warmUpUDP(WiFiUDP& udp, IPAddress serverIP, unsigned int serverPort) {
  udp.beginPacket(serverIP, serverPort);
  udp.write((uint8_t)0);  // ダミーデータ送信
  udp.endPacket();
  delay(50);  // 少し待機
}