#ifndef AGENT_CONFIG_H
#define AGENT_CONFIG_H

#include <WiFi.h>
#include <WiFiUdp.h>
#pragma once
#include <Arduino.h>

// サーバから omega, kappa, alpha を取得する関数
void requestParametersFromServer(WiFiUDP &udp, IPAddress serverIP, unsigned int serverPort, float &omega, float &kappa, float &alpha);

// agent_id をファイルから読み取る関数
int readAgentIdFromFile();

#endif  // AGENT_CONFIG_H