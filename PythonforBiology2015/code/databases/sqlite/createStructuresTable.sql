create table structures (
  id INTEGER,
  structureName TEXT NOT NULL,
  structurePdbId TEXT NOT NULL,
  structureModel INTEGER NOT NULL,
  chainMolType TEXT NOT NULL,
  chainCode TEXT NOT NULL,
  residueSeqId INTEGER NOT NULL,
  residueCode TEXT,
  atomName TEXT NOT NULL,
  atomX FLOAT NOT NULL,
  atomY FLOAT NOT NULL,
  atomZ FLOAT NOT NULL,
  atomElement TEXT NOT NULL,
  PRIMARY KEY (id)
);

