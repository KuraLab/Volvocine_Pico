#ifndef AGENT_CONFIG_H
#define AGENT_CONFIG_H

#include <WiFi.h>
#include <WiFiUdp.h>
#pragma once
#include <Arduino.h>

// サーバからパラメータを取得し、omega と kappa を更新する関数
// 引数:
//   udpParam         : パラメータ要求用のWiFiUDPインスタンス
//   serverIP         : サーバのIPアドレス
//   paramServerPort  : サーバ側パラメータ要求ポート番号
//   omega            : 取得した角速度（更新用、参照渡し）
//   kappa            : 取得したゲイン（更新用、参照渡し）
bool requestParameters(int agent_id, WiFiUDP &udpParam, IPAddress serverIP, unsigned int paramServerPort, float &omega, float &kappa);

// agent_id をファイルから読み取る関数
int readAgentIdFromFile();

#endif  // AGENT_CONFIG_H