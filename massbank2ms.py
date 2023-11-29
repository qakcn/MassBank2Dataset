from functions import *

conn = MySQLdb.connect(host='192.168.11.50', user='username', password='password', database='MassBank')
mspath = Path('./ms')

get_compounds(conn, mspath)
