#include "agent_config.h"
#include <LittleFS.h>

int readAgentIdFromFile() {
    // Pico(W) 版 LittleFS では引数なしで begin()
    if(!LittleFS.begin()){
        Serial.println("LittleFS Mount Failed!");
        return 0; // デフォルト値
    }

    File f = LittleFS.open("/config.txt", "r");
    if(!f){
        Serial.println("Failed to open /config.txt");
        return 0; // デフォルト値
    }

    // シングルクォートに '\n' と書く
    // ( '\\n' はマルチキャラクタリテラル扱いで警告が出る )
    String line = f.readStringUntil('\n');
    f.close();

    return line.toInt(); // ファイルの値をintに変換して返す
}

void requestParametersFromServer(WiFiUDP &udp, IPAddress serverIP, unsigned int serverPort, float &omega, float &kappa, float &alpha) {
  // リクエスト送信
  udp.beginPacket(serverIP, serverPort);
  udp.write("REQUEST_PARAMS");  // パラメータリクエスト用の識別文字列
  udp.endPacket();

  // 応答待機
  unsigned long startTime = millis();
  while (millis() - startTime < 2000) {  // 最大2秒待機
    int packetSize = udp.parsePacket();
    if (packetSize) {
      char buffer[128];
      int len = udp.read(buffer, sizeof(buffer) - 1);
      if (len > 0) {
        buffer[len] = '\0';  // 文字列終端を追加
        sscanf(buffer, "omega:%f,kappa:%f,alpha:%f", &omega, &kappa, &alpha);
        Serial.printf("[INFO] Received parameters: omega=%.2f, kappa=%.2f, alpha=%.2f\n", omega, kappa, alpha);
        return;
      }
    }
    delay(100);  // 少し待機
  }

  Serial.println("[WARN] Parameter request timed out.");
}
