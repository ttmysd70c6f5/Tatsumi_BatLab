const int pinDoor = 6;  //Pin Reed switch
const int pinOPEN    = 9;  //Pin LED for OPEN (red)
const int pinCLOSE    = 12;  //Pin LED for CLOSE (green)
int door = 0;
void setup()
{
  Serial.begin(9600); // start serial communication
  // setup pins
  pinMode(pinOPEN, OUTPUT);
  pinMode(pinCLOSE, OUTPUT);
  pinMode(pinDoor, INPUT);
}
void loop()
{
  // Update door LEDs
  door = digitalRead(pinDoor);  //Read pin for reed switch
  if (door == HIGH)
  {
    digitalWrite(pinOPEN, HIGH);
    digitalWrite(pinCLOSE, LOW);
  }
  else
  {
    digitalWrite(pinOPEN, LOW);
    digitalWrite(pinCLOSE, HIGH);
  }

  // Send door status and TTL signals to PC
  int ttl = analogRead(0);  // Read analog pin for TTLs
  door = digitalRead(pinDoor);  //Read pin for reed switch
  unsigned long time = millis();  // Time passed since starting arduino

  char buff[18];
  sprintf(buff, "%lu,%i,%i", time, ttl, door);
  Serial.println(buff);
  delay(1000);
}
