/*
 * adaptive_control.ino
 * 
 * EECS 127 Spring 2020
 * Vade Shah
 * 
 * Based heavily upon code developed by John Maidens, Emily Naviasky, Nathaniel Mailoa, and Andrew Blatner
 * from UC Berkeley's EE 16B
 * 
 */

#include "MatrixMath.h"

#define LEFT_MOTOR                  P2_0
#define LEFT_ENCODER                P6_1
#define RIGHT_MOTOR                 P1_5
#define RIGHT_ENCODER               P6_2

#define RUN_TIME                    (5*1000)                     // The car is only run for 50-100ms intervals
#define SAMPLING_INTERVAL           100                          // totaling 5 seconds. See the report for more
#define SAMPLE_LEN                  (RUN_TIME/SAMPLING_INTERVAL) // information.

#define JOLT_STEPS                  2

MatrixMath m;

int step_num = 0;
volatile boolean do_loop = 0; // timer signal to increment timestep

int16_t deltaArr[SAMPLE_LEN] = {0};
uint16_t lpos[SAMPLE_LEN] = {0};
uint16_t rpos[SAMPLE_LEN] = {0};
uint8_t lpwm[SAMPLE_LEN] = {0};
uint8_t rpwm[SAMPLE_LEN] = {0};

typedef struct encoder {
  int pin;
  int pos;
  bool level;
  int avg;
} encoder_t;

encoder_t left_encoder = {LEFT_ENCODER, 0, LOW, 0};
encoder_t right_encoder = {RIGHT_ENCODER, 0, LOW, 0};

// these initial values were obtained via least squares previously
float theta_left = 0.06899;
float theta_right = 0.2982;
float beta_left = -51.46;
float beta_right = -19.55;
float v_star = 64.8;

// PWM inputs to jolt the car straight
int left_jolt = 235;
int right_jolt= 240;

// Control gains
float k_left = 0.5;
float k_right = 0.5;

float delta_ss = 2; // the steady states difference between the left and right wheels

void setup(void) {
  Serial.begin(38400);

  pinMode(LEFT_MOTOR, OUTPUT);
  pinMode(LEFT_ENCODER, INPUT);
  pinMode(RIGHT_MOTOR, OUTPUT);
  pinMode(RIGHT_ENCODER, INPUT);
  pinMode(RED_LED, OUTPUT);
  pinMode(GREEN_LED, OUTPUT);

  write_pwm(0, 0); // Turn off motors
  delay(2000); // Wait 2 seconds to put down car
  reset_blinker(); // Blink lights to indicate car is running
  setTimer(); // Set timer for timestep
}

void loop(void) {
  // put your main code here, to run repeatedly: 

  //alpha = *ptrResult;
  //beta = *(ptrResult + 1);

  check_encoders();
  if (do_loop) {
    // Apply maximum input for a short time to start motors
    if (step_num < JOLT_STEPS) {
      write_pwm(left_jolt, right_jolt);
      step_num++;
    }
    // If not done running
    else if (step_num == 2 or step_num == 3 or step_num == 4) {
      // Save positions because _left_position and _right_position
      // can change in the middle of one loop.
      int left_position = left_encoder.pos;
      int right_position = right_encoder.pos;

      float delta = left_position - right_position + delta_ss;

      // Drive straight using feedback
      // Compute the needed pwm values for each wheel using delta and v_star
      int left_cur_pwm = driveStraight_left(v_star, delta);
      int right_cur_pwm = driveStraight_right(v_star, delta);
      write_pwm(left_cur_pwm, right_cur_pwm);

      lpos[step_num] = left_position;
      rpos[step_num] = right_position;
      deltaArr[step_num] = delta;
      lpwm[step_num] = left_cur_pwm;
      rpwm[step_num] = right_cur_pwm;

      step_num++;
    }
    else if (step_num < SAMPLE_LEN) {

      // Save positions because _left_position and _right_position
      // can change in the middle of one loop.
      int left_position = left_encoder.pos;
      int right_position = right_encoder.pos;

      float delta = left_position - right_position + delta_ss;

      // Drive straight using feedback
      // Compute the needed pwm values for each wheel using delta and v_star
      // Compute theta_right, beta_right, theta_left, and beta_left using least squares
      
      int h = step_num - 3;

      // loop computes A
      float AL[2 * h];
      for (int i = 3; i < step_num; i++) {
        AL[2 * i + 0 - 6] = lpwm[i - 1];
        AL[2 * i + 1 - 6] = -1;
      }

      // loop computes B
      float bL[h];
      for (int i = 3; i < step_num; i++) {
       bL[i - 3] = lpos[i] - lpos[i - 1];
      }
  
      float DL[2 * h]; // At
      float EL[4]; // At*A and (At*A)-1
      float FL[2 * h]; // (At*A)-1 * At
      float GL[2];
      m.MatrixMath::MatrixCopy(AL, h, 2, DL);
      m.MatrixMath::MatrixTranspose(AL, h, 2, DL);
      m.MatrixMath::MatrixMult(DL, AL, 2, h, 2, EL);
      m.MatrixMath::MatrixInvert(EL, 2);
      m.MatrixMath::MatrixMult(EL, DL, 2, 2, h, FL);
      m.MatrixMath::MatrixMult(FL, bL, 2, h, 1, GL);
      theta_left = GL[0];
      beta_left = GL[1];

      //calculations for right parameters

      // loop computes A
      float AR[2 * h];
      for (int i = 3; i < step_num; i++) {
        AR[2 * i + 0 - 6] = rpwm[i - 1];
        AR[2 * i + 1 - 6] = -1;
      }

      // loop computes B
      float bR[h];
      for (int i = 3; i < step_num; i++) {
        bR[i - 3] = rpos[i] - rpos[i - 1];
      }

      float DR[2 * h]; // At
      float ER[4]; // At*A and (At*A)-1
      float FR[2 * h]; // (At*A)-1 * At
      float GR[2];
      m.MatrixMath::MatrixCopy(AR, h, 2, DR);
      m.MatrixMath::MatrixTranspose(AR, h, 2, DR);
      m.MatrixMath::MatrixMult(DR, AR, 2, h, 2, ER);
      m.MatrixMath::MatrixInvert(ER, 2);
      m.MatrixMath::MatrixMult(ER, DR, 2, 2, h, FR);
      m.MatrixMath::MatrixMult(FR, bL, 2, h, 1, GR);
      theta_right = GR[0];
      beta_right = GR[1];
      
      int left_cur_pwm = driveStraight_left_ac(v_star, delta, beta_left, theta_left);
      int right_cur_pwm = driveStraight_right_ac(v_star, delta, beta_right, theta_right);
      write_pwm(left_cur_pwm, right_cur_pwm);

      lpos[step_num] = left_position;
      rpos[step_num] = right_position;
      deltaArr[step_num] = delta;
      lpwm[step_num] = left_cur_pwm;
      rpwm[step_num] = right_cur_pwm;

      step_num++;
    }

    else { // When step_num has reached SAMPLE_LEN
      // Turn off motors
      write_pwm(0, 0);

      // Print out result
      Serial.println("Start");
      Serial.println("delta - left pos - right pos - left pwm - right pwm");
      for (int i = 0; i < SAMPLE_LEN; i++) {
        Serial.print(deltaArr[i]);
        Serial.print(',');
        Serial.print(lpos[i]);
        Serial.print(',');
        Serial.print(rpos[i]);
        Serial.print(',');
        Serial.print(lpwm[i]);
        Serial.print(',');
        Serial.print(rpwm[i]);
        Serial.print('\n');
      }
    }
    do_loop = 0;
  }
}

/*----------------------------------*/
/*     Custom helper functions      */
/*----------------------------------*/

float* lstsq(uint8_t u[], uint16_t p[], int height){
  // given a least squares setup Ax = b, this function computes (At*A)-1*At*b
  int h = 2 * height - 6;
  float* A = mat(u, height); // A
  float* b = vec(p, height); // b
  float D[h]; // At
  float E[4]; // At*A and (At*A)-1
  float F[h]; // (At*A)-1 * At
  float G[2];
  float* ptrD = D;
  float* ptrE = E;
  float* ptrF = F;
  float* ptrG = G;

  m.MatrixMath::MatrixCopy(A, h, 2, ptrD);
  m.MatrixMath::MatrixTranspose(A, h, 2, ptrD);
  m.MatrixMath::MatrixMult(ptrD, A, 2, h, 2, ptrE);
  m.MatrixMath::MatrixInvert(ptrE, 2);
  m.MatrixMath::MatrixMult(ptrE, ptrD, 2, 2, h, ptrF);
  m.MatrixMath::MatrixMult(ptrF, b, 2, h, 1, ptrG);
  return ptrG;
}

float driveStraight_left_ac(float v_star, float delta, float beta_left, float theta_left) {
  return (v_star + beta_left) / theta_left - k_left * delta / theta_left;
}

float driveStraight_right_ac(float v_star, float delta, float beta_right, float theta_right) {
  return (v_star + beta_right) / theta_right + k_right * delta / theta_right;
}

float* mat(uint8_t u[], int height) {
  float m[2 * height - 6];
  for (int i = 3; i < height; i++) {
    m[2 * i + 0 - 6] = u[i - 1];
    m[2 * i + 1 - 6] = -1;
  }
  float* ptr = m;
  return ptr;
}

float* vec(uint16_t p[], int height) {
  float v[height - 3];
  for (int i = 3; i < height; i++) {
    v[i - 3] = p[i] - p[i - 1];
  }
  float* ptr = v;
  return ptr;
}

/*---------------------------*/
/*     Helper functions      */
/*---------------------------*/

void write_pwm(int pwm_left, int pwm_right) {
  analogWrite(LEFT_MOTOR, (int) min(max(0, pwm_left), 255));
  analogWrite(RIGHT_MOTOR, (int) min(max(0, pwm_right), 255));
}

void reset_blinker(void) {
  digitalWrite(RED_LED, HIGH);
  delay(100);
  digitalWrite(RED_LED, LOW);
  digitalWrite(GREEN_LED, HIGH);
  delay(100);
  digitalWrite(RED_LED, HIGH);
  digitalWrite(GREEN_LED, LOW);
  delay(100);
  digitalWrite(RED_LED, LOW);
  digitalWrite(GREEN_LED, HIGH);
  delay(100);
  digitalWrite(GREEN_LED, LOW);
}

float driveStraight_left(float v_star, float delta) {
  return (v_star + beta_left) / theta_left - k_left * delta / theta_left;
}

float driveStraight_right(float v_star, float delta) {
  return (v_star + beta_right) / theta_right + k_right * delta / theta_right;
}

/*---------------------------*/
/*    Interrupt functions    */
/*---------------------------*/

#define AVG_DECAY_RATE              0.3
#define LOW_THRESH                  ((int) (0.2*4096))
#define HIGH_THRESH                 ((int) (0.8*4096))

void check_encoder(encoder_t* enc) {
  int new_val = analogRead(enc->pin);
  enc->avg = (int) (AVG_DECAY_RATE*enc->avg + (1 - AVG_DECAY_RATE)*new_val);
  if ((enc->level == LOW && HIGH_THRESH < enc->avg) ||
      (enc->level == HIGH && enc->avg < LOW_THRESH)) {
    enc->pos++;
    enc->level = !enc->level;
  }
}

void check_encoders(void) {
  check_encoder(&left_encoder);
  check_encoder(&right_encoder);
}

// Set timer for timestep; use A2 since A0 & A1 are used by PWM
void setTimer(void) {
  TA2CCR0 = (unsigned int) (32.768*SAMPLING_INTERVAL); // set the timer based on 32kHz clock
  TA2CCTL0 = CCIE; // enable interrupts for Timer A
  __bis_SR_register(GIE);
  TA2CTL = TASSEL_1 + MC_1 + TACLR + ID_0;
}

// ISR for timestep
#pragma vector=TIMER2_A0_VECTOR    // Timer A ISR
__interrupt void Timer2_A0_ISR(void) {
  do_loop = 1;
}
