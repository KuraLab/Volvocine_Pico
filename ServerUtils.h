#ifndef SERVER_UTILS_H
#define SERVER_UTILS_H

#include <WiFiUdp.h>

// サーバーの応答を確認する関数
bool isServerReady(WiFiUDP& udp, IPAddress serverIP, unsigned int serverPort);

// UDP通信のウォームアップ関数
void warmUpUDP(WiFiUDP& udp, IPAddress serverIP, unsigned int serverPort);

#endif // SERVER_UTILS_H