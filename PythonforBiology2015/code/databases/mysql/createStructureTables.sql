create table structure (
  id INTEGER AUTO_INCREMENT,
  name TEXT NOT NULL,
  pdbId TEXT NOT NULL,
  conformation INTEGER NOT NULL,
  PRIMARY KEY (id)
);

create table chain (
  id INTEGER AUTO_INCREMENT,
  structureId INTEGER NOT NULL,
  molType TEXT NOT NULL,
  code TEXT NOT NULL,
  PRIMARY KEY (id),
  FOREIGN KEY (structureId) REFERENCES structure(id)
);

create table residue (
  id INTEGER AUTO_INCREMENT,
  chainId INTEGER NOT NULL,
  seqId INTEGER NOT NULL,
  code TEXT,
  PRIMARY KEY (id),
  FOREIGN KEY (chainId) REFERENCES chain(id)
);

create table atom (
  id INTEGER AUTO_INCREMENT,
  residueId INTEGER NOT NULL,
  name TEXT NOT NULL,
  x FLOAT NOT NULL,
  y FLOAT NOT NULL,
  z FLOAT NOT NULL,
  element TEXT NOT NULL,
  PRIMARY KEY (id),
  FOREIGN KEY (residueId) REFERENCES residue(id)
);

