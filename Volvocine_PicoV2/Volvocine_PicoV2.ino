#include <Servo.h>
#include <Wire.h>
#include <math.h>
#include <stdio.h>

Servo myServo;

// Volvocine_Pico.ino と同じサーボピン
const int servoPin = 1;
const int switchPin = 4;
const int analogPin1 = 27;
const int analogPin2 = 28;

// 動作パラメータ
const int centerAngle = 90;
const int amplitude = 60;
const float frequencyHz = 1.25f;  // 0.5Hz = 2秒で1往復

// INA226設定
uint8_t ina226Addr = 0x40;
const float SHUNT_OHMS = 0.056f;  // 実測系のシャント抵抗値
unsigned long lastInaPrintMs = 0;
const unsigned long inaPrintIntervalMs = 2;
const uint8_t inaSamplesPerSend = 1;
const unsigned long inaReadErrorLogIntervalMs = 1000;

const bool enableLoopProfiler = false;
const unsigned long profilerReportIntervalMs = 5000;

float currentSummA = 0.0f;
float busSumV = 0.0f;
uint8_t inaSampleCount = 0;
unsigned long long timeSumUs = 0;

unsigned long long profSwitchUsSum = 0;
unsigned long long profServoUsSum = 0;
unsigned long long profInaUsSum = 0;
unsigned long long profInaReadUsSum = 0;
unsigned long long profTxUsSum = 0;
unsigned long long profLoopUsSum = 0;
unsigned long profSwitchUsMax = 0;
unsigned long profServoUsMax = 0;
unsigned long profInaUsMax = 0;
unsigned long profInaReadUsMax = 0;
unsigned long profTxUsMax = 0;
unsigned long profLoopUsMax = 0;
unsigned long profLoopCount = 0;
unsigned long lastProfReportMs = 0;

bool motionMode = false;
bool lastSwitchState = false;
unsigned long lastToggleTime = 0;
const unsigned long debounceMs = 200;
float phase = 0.0f;
unsigned long lastPhaseUpdateMs = 0;
bool inaReady = false;
unsigned long lastInaReadErrorLogMs = 0;

bool inaWriteReg16(uint8_t reg, uint16_t value) {
  Wire1.beginTransmission(ina226Addr);
  Wire1.write(reg);
  Wire1.write((uint8_t)(value >> 8));
  Wire1.write((uint8_t)(value & 0xFF));
  return Wire1.endTransmission() == 0;
}

bool inaReadReg16(uint8_t reg, uint16_t &value) {
  Wire1.beginTransmission(ina226Addr);
  Wire1.write(reg);
  if (Wire1.endTransmission(false) != 0) {
    return false;
  }

  if (Wire1.requestFrom((int)ina226Addr, 2) != 2) {
    return false;
  }

  value = ((uint16_t)Wire1.read() << 8) | (uint16_t)Wire1.read();
  return true;
}

void initIna226() {
  // 高速設定: 平均1回、変換時間140us、連続シャント/バス変換
  // 0x0007 = AVG=1, VBUSCT=140us, VSHCT=140us, MODE=111
  inaWriteReg16(0x00, 0x0007);
}

bool detectIna226Address() {
  if (Serial) Serial.println("[INA226] scanning I2C1 addresses 0x40-0x4F...");
  for (uint8_t addr = 0x40; addr <= 0x4F; addr++) {
    Wire1.beginTransmission(addr);
    uint8_t err = Wire1.endTransmission();
    if (err == 0) {
      ina226Addr = addr;
      if (Serial) Serial.print("[INA226] device found at 0x");
      if (Serial) Serial.println(ina226Addr, HEX);
      return true;
    }
  }
  if (Serial) Serial.println("[INA226] no device found in 0x40-0x4F");
  return false;
}

bool readIna226Measurement(float &currentmA, float &busVoltV) {
  uint16_t shuntRawU16 = 0;
  uint16_t busRawU16 = 0;

  if (!inaReadReg16(0x01, shuntRawU16) || !inaReadReg16(0x02, busRawU16)) {
    unsigned long nowMs = millis();
    if (nowMs - lastInaReadErrorLogMs >= inaReadErrorLogIntervalMs) {
      lastInaReadErrorLogMs = nowMs;
        if (Serial) Serial.println("[INA226] read error (check wiring/pull-up/address)");
    }
    return false;
  }

  int16_t shuntRaw = (int16_t)shuntRawU16;
  float shuntVoltV = (float)shuntRaw * 2.5e-6f;     // 2.5uV/bit
  float currentA = shuntVoltV / SHUNT_OHMS;
  currentmA = currentA * 1000.0f;
  busVoltV = (float)busRawU16 * 1.25e-3f;     // 1.25mV/bit

  return true;
}

bool sendIna226Measurement(float currentmA, float busVoltV, uint32_t timestampUs, uint16_t adc1, uint16_t adc2) {
  // 送信でループが止まらないよう、空きが足りない時はこのサンプルをスキップする
  char line[80];
  int lineLen = snprintf(line, sizeof(line), "%.2f,%.3f,%lu,%u,%u\n",
                         currentmA, busVoltV, (unsigned long)timestampUs,
                         (unsigned int)adc1, (unsigned int)adc2);
  if (lineLen <= 0 || lineLen >= (int)sizeof(line)) {
    return false;
  }

  if (Serial.availableForWrite() < lineLen) {
    return false;
  }

  size_t written = Serial.write((const uint8_t *)line, (size_t)lineLen);
  return written == (size_t)lineLen;
}

void setup() {
  Serial.begin(921600);
  delay(200);

  // I2C1をGP6(SDA), GP7(SCL)に明示的に割り当て
  Wire1.setSDA(6);
  Wire1.setSCL(7);
  Wire1.begin();
  Wire1.setClock(400000);
  inaReady = detectIna226Address();
  if (inaReady) {
    initIna226();
  }

  pinMode(switchPin, INPUT);
  analogReadResolution(12);
  myServo.attach(servoPin);
  myServo.write(centerAngle);
  lastPhaseUpdateMs = millis();
  lastProfReportMs = millis();

  if (Serial) Serial.println("[INFO] Servo + INA226 monitor start");
}

void loop() {
  unsigned long loopStartUs = micros();

  unsigned long sectionStartUs = micros();
  bool currentSwitchState = digitalRead(switchPin);
  if (currentSwitchState && !lastSwitchState) {
    unsigned long now = millis();
    if (now - lastToggleTime >= debounceMs) {
      motionMode = !motionMode;
      lastToggleTime = now;

      if (motionMode) {
        // 再開時は位相を維持したまま、時間差分だけリセットしてジャンプを防ぐ
        lastPhaseUpdateMs = now;
      }
    }
  }
  lastSwitchState = currentSwitchState;
  unsigned long switchUs = micros() - sectionStartUs;

  sectionStartUs = micros();
  if (motionMode) {
    unsigned long now = millis();
    float dt = (now - lastPhaseUpdateMs) / 1000.0f;
    lastPhaseUpdateMs = now;

    phase += 2.0f * PI * frequencyHz * dt;
    if (phase >= 2.0f * PI) {
      phase = fmodf(phase, 2.0f * PI);
    }

    float wave = sinf(phase);
    int angle = (int)(centerAngle + amplitude * wave);

    myServo.write(angle);
  }
  unsigned long servoUs = micros() - sectionStartUs;

  sectionStartUs = micros();
  unsigned long inaReadUs = 0;
  unsigned long txUs = 0;
  unsigned long nowMs = millis();
  if (nowMs - lastInaPrintMs >= inaPrintIntervalMs) {
    lastInaPrintMs = nowMs;
    if (inaReady) {
      float currentmA = 0.0f;
      float busVoltV = 0.0f;

      unsigned long readStartUs = micros();
      if (readIna226Measurement(currentmA, busVoltV)) {
        inaReadUs = micros() - readStartUs;
        uint16_t adc1 = (uint16_t)analogRead(analogPin1);
        uint16_t adc2 = (uint16_t)analogRead(analogPin2);
        currentSummA += currentmA;
        busSumV += busVoltV;
        timeSumUs += readStartUs;
        inaSampleCount++;

        if (inaSampleCount >= inaSamplesPerSend) {
          float invN = 1.0f / (float)inaSampleCount;
          uint32_t avgTimestampUs = (uint32_t)(timeSumUs / (unsigned long long)inaSampleCount);
          unsigned long txStartUs = micros();
          sendIna226Measurement(currentSummA * invN, busSumV * invN, avgTimestampUs, adc1, adc2);
          txUs = micros() - txStartUs;
          // 送信成否に関わらず蓄積をクリアし、古いデータの滞留を防ぐ
          currentSummA = 0.0f;
          busSumV = 0.0f;
          timeSumUs = 0;
          inaSampleCount = 0;
        }
      } else {
        inaReadUs = micros() - readStartUs;
      }
    }
  }
  unsigned long inaUs = micros() - sectionStartUs;

  unsigned long loopUs = micros() - loopStartUs;

  if (enableLoopProfiler) {
    profLoopCount++;
    profSwitchUsSum += switchUs;
    profServoUsSum += servoUs;
    profInaUsSum += inaUs;
    profInaReadUsSum += inaReadUs;
    profTxUsSum += txUs;
    profLoopUsSum += loopUs;

    if (switchUs > profSwitchUsMax) profSwitchUsMax = switchUs;
    if (servoUs > profServoUsMax) profServoUsMax = servoUs;
    if (inaUs > profInaUsMax) profInaUsMax = inaUs;
    if (inaReadUs > profInaReadUsMax) profInaReadUsMax = inaReadUs;
    if (txUs > profTxUsMax) profTxUsMax = txUs;
    if (loopUs > profLoopUsMax) profLoopUsMax = loopUs;

    unsigned long nowMsForProf = millis();
    if (nowMsForProf - lastProfReportMs >= profilerReportIntervalMs && profLoopCount > 0) {
      lastProfReportMs = nowMsForProf;
      float avgLoopUs = (float)profLoopUsSum / (float)profLoopCount;
      float controlPeriodMs = avgLoopUs / 1000.0f;
      float controlHz = avgLoopUs > 0.0f ? 1000000.0f / avgLoopUs : 0.0f;

      if (Serial) Serial.print("[PROF] avg_us sw=");
      if (Serial) Serial.print((float)profSwitchUsSum / (float)profLoopCount, 2);
      if (Serial) Serial.print(", servo=");
      if (Serial) Serial.print((float)profServoUsSum / (float)profLoopCount, 2);
      if (Serial) Serial.print(", ina=");
      if (Serial) Serial.print((float)profInaUsSum / (float)profLoopCount, 2);
      if (Serial) Serial.print(", ina_read=");
      if (Serial) Serial.print((float)profInaReadUsSum / (float)profLoopCount, 2);
      if (Serial) Serial.print(", tx=");
      if (Serial) Serial.print((float)profTxUsSum / (float)profLoopCount, 2);
      if (Serial) Serial.print(", loop=");
      if (Serial) Serial.print(avgLoopUs, 2);
      if (Serial) Serial.print(" | ctrl=");
      if (Serial) Serial.print(controlPeriodMs, 4);
      if (Serial) Serial.print(" ms (");
      if (Serial) Serial.print(controlHz, 1);
      if (Serial) Serial.print(" Hz)");
      if (Serial) Serial.print(" | max_us sw=");
      if (Serial) Serial.print(profSwitchUsMax);
      if (Serial) Serial.print(", servo=");
      if (Serial) Serial.print(profServoUsMax);
      if (Serial) Serial.print(", ina=");
      if (Serial) Serial.print(profInaUsMax);
      if (Serial) Serial.print(", ina_read=");
      if (Serial) Serial.print(profInaReadUsMax);
      if (Serial) Serial.print(", tx=");
      if (Serial) Serial.print(profTxUsMax);
      if (Serial) Serial.print(", loop=");
      if (Serial) Serial.println(profLoopUsMax);

      profSwitchUsSum = 0;
      profServoUsSum = 0;
      profInaUsSum = 0;
      profInaReadUsSum = 0;
      profTxUsSum = 0;
      profLoopUsSum = 0;
      profSwitchUsMax = 0;
      profServoUsMax = 0;
      profInaUsMax = 0;
      profInaReadUsMax = 0;
      profTxUsMax = 0;
      profLoopUsMax = 0;
      profLoopCount = 0;
    }
  }

  delay(0); // 可能な限り高速にループ
}
