CREATE TABLE Person (
  firstName VARCHAR(30) NOT NULL,
  lastName VARCHAR(30) NOT NULL,
  birthYear INT,
  PRIMARY KEY (firstName, lastName)
);

