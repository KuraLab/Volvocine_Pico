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

void requestParametersFromServer(WiFiUDP &udp, IPAddress serverIP, unsigned int serverPort, int agent_id, float &omega, float &kappa, float &alpha, float &servoCenter, float &servoAmplitude, int &stopAgentId, int &stopDelaySeconds) {
  // デバッグ情報を含むリクエスト文字列を作成
  int analogValue26 = analogRead(26);  // 26ピンのアナログ入力値を取得
  char requestBuffer[128]; // バッファサイズを拡張して新しいパラメータに対応
  snprintf(requestBuffer, sizeof(requestBuffer), "REQUEST_PARAMS,id:%d,analog26:%d", agent_id, analogValue26);

  // リクエスト送信
  udp.beginPacket(serverIP, serverPort);
  udp.write(requestBuffer);  // デバッグ情報付きリクエスト
  udp.endPacket();

  // 応答待機
  unsigned long startTime = millis();
  while (millis() - startTime < 2000) {  // 最大2秒待機
    int packetSize = udp.parsePacket();
    if (packetSize) {
      char buffer[256]; // 応答バッファサイズを拡張
      int len = udp.read(buffer, sizeof(buffer) - 1);
      if (len > 0) {
        buffer[len] = '\0';  // 文字列終端を追加
        // sscanfでのパースを修正。新しいパラメータに対応。
        // 例: "omega:3.14,kappa:10.0,alpha:-1.57,center:110.0,amplitude:60.0,stop_id:2,stop_delay:100"
        int parsed_count = sscanf(buffer, "omega:%f,kappa:%f,alpha:%f,center:%f,amplitude:%f,stop_id:%d,stop_delay:%d", &omega, &kappa, &alpha, &servoCenter, &servoAmplitude, &stopAgentId, &stopDelaySeconds);
        if (parsed_count == 7) { // 7つのパラメータすべてがパースできたか確認
            Serial.printf("[INFO] Received parameters: omega=%.2f, kappa=%.2f, alpha=%.2f, center=%.1f, amplitude=%.1f, stop_id=%d, stop_delay=%d\n", omega, kappa, alpha, servoCenter, servoAmplitude, stopAgentId, stopDelaySeconds);
        } else {
            Serial.printf("[WARN] Failed to parse all parameters. Received: %s (parsed %d)\n", buffer, parsed_count);
        }
        return;
      }
    }
    delay(100);  // 少し待機
  }

  Serial.println("[WARN] Parameter request timed out.");
}
