#include "agent_config.h"
#include <LittleFS.h>
#include <math.h>
#include <string.h>

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

void requestParametersFromServer(WiFiUDP &udp, IPAddress serverIP, unsigned int serverPort, int agent_id, float &omega, float &kappa, float &alpha, float &servoCenter, float &servoAmplitude, int &stopAgentId, int &stopDelaySeconds, int &prcHarmonics, float *prcCosCoeffs, float *prcSinCoeffs, int prcMaxHarmonics) {
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
      char buffer[1024]; // フーリエ係数を含む応答を受けるため拡張
      int len = udp.read(buffer, sizeof(buffer) - 1);
      if (len > 0) {
        buffer[len] = '\0';  // 文字列終端を追加
        char originalBuffer[1024];
        strncpy(originalBuffer, buffer, sizeof(originalBuffer) - 1);
        originalBuffer[sizeof(originalBuffer) - 1] = '\0';

        int parsedBaseFields = 0;
        int parsedPrcFields = 0;
        bool gotPrcN = false;
        int receivedPrcN = 0;
        float tempCos[16] = {0.0f};
        float tempSin[16] = {0.0f};

        char *token = strtok(buffer, ",");
        while (token != nullptr) {
          float fval = 0.0f;
          int ival = 0;
          int idx = 0;

          if (sscanf(token, "omega:%f", &fval) == 1) {
            omega = fval;
            parsedBaseFields++;
          } else if (sscanf(token, "kappa:%f", &fval) == 1) {
            kappa = fval;
            parsedBaseFields++;
          } else if (sscanf(token, "alpha:%f", &fval) == 1) {
            alpha = fval;
            parsedBaseFields++;
          } else if (sscanf(token, "center:%f", &fval) == 1) {
            servoCenter = fval;
            parsedBaseFields++;
          } else if (sscanf(token, "amplitude:%f", &fval) == 1) {
            servoAmplitude = fval;
            parsedBaseFields++;
          } else if (sscanf(token, "stop_id:%d", &ival) == 1) {
            stopAgentId = ival;
            parsedBaseFields++;
          } else if (sscanf(token, "stop_delay:%d", &ival) == 1) {
            stopDelaySeconds = ival;
            parsedBaseFields++;
          } else if (sscanf(token, "prc_n:%d", &ival) == 1) {
            receivedPrcN = ival;
            gotPrcN = true;
            parsedPrcFields++;
          } else if (sscanf(token, "prc_a%d:%f", &idx, &fval) == 2) {
            if (idx >= 1 && idx <= prcMaxHarmonics && idx < (int)(sizeof(tempCos) / sizeof(tempCos[0]))) {
              tempCos[idx] = fval;
              parsedPrcFields++;
            }
          } else if (sscanf(token, "prc_b%d:%f", &idx, &fval) == 2) {
            if (idx >= 1 && idx <= prcMaxHarmonics && idx < (int)(sizeof(tempSin) / sizeof(tempSin[0]))) {
              tempSin[idx] = fval;
              parsedPrcFields++;
            }
          }

          token = strtok(nullptr, ",");
        }

        // PRCが来ない場合は従来式 cos(psi-alpha) と等価な1次フーリエ係数へフォールバック。
        for (int n = 0; n <= prcMaxHarmonics; n++) {
          prcCosCoeffs[n] = 0.0f;
          prcSinCoeffs[n] = 0.0f;
        }

        if (gotPrcN && receivedPrcN > 0 && parsedPrcFields > 1) {
          if (receivedPrcN > prcMaxHarmonics) {
            receivedPrcN = prcMaxHarmonics;
          }
          prcHarmonics = receivedPrcN;
          for (int n = 1; n <= prcHarmonics; n++) {
            prcCosCoeffs[n] = tempCos[n];
            prcSinCoeffs[n] = tempSin[n];
          }
        } else {
          prcHarmonics = 1;
          prcCosCoeffs[1] = cosf(alpha);
          prcSinCoeffs[1] = sinf(alpha);
        }

        if (parsedBaseFields >= 7) {
            Serial.printf("[INFO] Received parameters: omega=%.2f, kappa=%.2f, alpha=%.2f, center=%.1f, amplitude=%.1f, stop_id=%d, stop_delay=%d, prc_n=%d (parsed_prc=%d)\n", omega, kappa, alpha, servoCenter, servoAmplitude, stopAgentId, stopDelaySeconds, prcHarmonics, parsedPrcFields);
        } else {
          Serial.printf("[WARN] Failed to parse all base parameters. Received: %s (parsed_base=%d, parsed_prc=%d)\n", originalBuffer, parsedBaseFields, parsedPrcFields);
        }
        return;
      }
    }
    delay(100);  // 少し待機
  }

  Serial.println("[WARN] Parameter request timed out. Keeping last received PRC.");
}
