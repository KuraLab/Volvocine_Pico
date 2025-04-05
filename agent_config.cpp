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

bool requestParameters(int agent_id, WiFiUDP &udpParam, IPAddress serverIP, unsigned int paramServerPort, float &omega, float &kappa) {
  // エージェント固有のパラメータ要求メッセージを作成
  char reqMsg[50];
  sprintf(reqMsg, "REQUEST_PARAM,agent=%d", agent_id);
  
  // サーバに要求メッセージを送信
  udpParam.beginPacket(serverIP, paramServerPort);
  udpParam.print(reqMsg);
  udpParam.endPacket();

  // 応答を待つ（タイムアウト2000ms）
  unsigned long startTime = millis();
  while (millis() - startTime < 2000) {
    int packetSize = udpParam.parsePacket();
    if (packetSize) {
      char buff[100];
      int len = udpParam.read(buff, sizeof(buff) - 1);
      if (len > 0) {
        buff[len] = '\0';
      }
      String reply = String(buff);
      // 例: "PARAM,omega=0.1,kappa=1.5"
      if (reply.startsWith("PARAM,")) {
        int idxOmega = reply.indexOf("omega=");
        int idxKappa = reply.indexOf("kappa=");
        if (idxOmega != -1 && idxKappa != -1) {
          int commaAfterOmega = reply.indexOf(',', idxOmega);
          String omegaStr;
          if (commaAfterOmega == -1) {
            omegaStr = reply.substring(idxOmega + 6);
          } else {
            omegaStr = reply.substring(idxOmega + 6, commaAfterOmega);
          }
          String kappaStr = reply.substring(idxKappa + 6);
          omega = omegaStr.toFloat();
          kappa = kappaStr.toFloat();
          Serial.println("Received parameters from server:");
          Serial.print("  omega = ");
          Serial.println(omega, 6);
          Serial.print("  kappa = ");
          Serial.println(kappa, 6);
          return true;
        }
      }
    }
    delay(10);
  }
  Serial.println("Parameter request timed out.");
  return false;
}
