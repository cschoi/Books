def formatStatement(text, placeHolder):

  if placeHolder == '%s':
    return text

  numInserts = text.count('%s')

  return text % (numInserts*(placeHolder,))

def countChainMolType(connection, molType, placeHolder='%s'):

  cursor = connection.cursor()

  try:
    # get matching chain records from database

    stmt = "select * from chain where molType=%s"
    stmt = formatStatement(stmt, placeHolder)

    cursor.execute(stmt, (molType,))
    result = cursor.fetchall()

    count = len(result)

  finally:
    cursor.close()

  return count

if __name__ == '__main__':

  import sys

  if len(sys.argv) != 3:
    print('need to specify database, molType')
    sys.exit(1)

  database = sys.argv[1]
  molType = sys.argv[2]

  import sqlite3
  connection = sqlite3.connect(database)
  placeHolder = '?'

  try:
    count = countChainMolType(connection, molType, placeHolder)
    print('Found %d chain(s)' % count)
  finally:
    connection.close()
