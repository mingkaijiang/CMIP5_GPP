# This is a simple, very basic example on how to connect to the sqlite cmip5_raijin_latest.db database and read data from it using the sqlite3 module
# the databse only table 'cmip5' columns are:
# id text, variable text, mip text, model text, experiment text, ensemble text, version text
# the id is actually the file path in unofficial-ESG-replica/tmp/tree

import sqlite3

# this opens a connection to the database
# NB if the database doesn't exists then it will create one
conn = sqlite3.connect('/g/data1/ua6/unofficial-ESG-replica/tmp/tree/cmip5_raijin_latest.db')

# From the connection, we get the cursor object. The cursor is used to traverse the records from the result set. We call the execute() method of the cursor and execute the SQL statement.
c = conn.cursor()
cursor = c.execute(''' SELECT * FROM cmip5 where variable==\"tas\" and experiment==\"rcp45\" and mip==\"Amon\" ''')
rows=cursor.fetchall()
# this will return the directory path of all matching ensembles
for row in rows:
    print row[0]

conn.close()
