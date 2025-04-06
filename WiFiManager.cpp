#include "WiFiManager.h"
#include <Arduino.h>

void connectToWiFi(const char* ssid, const char* password) {
  Serial.print("Connecting to WiFi...");
  WiFi.begin(ssid, password);

  // WiFi接続待機ループ
  while (WiFi.status() != WL_CONNECTED) {
    delay(500);
    Serial.print(".");
  }

  Serial.println("\nWiFi connected.");
  Serial.print("IP address: ");
  Serial.println(WiFi.localIP());
}