import sys
import re

class PaftolDpDetails(object):
    
    reDbusername = re.compile('username: *([^ ]+)')
    reDbpassword = re.compile('password: *([^ ]+)')
    reDbhost = re.compile('host: *([^ ]+)')
    reDbname = re.compile('dbname: *([^ ]+)')

    def __init__(self, detailsFile=None):
        if detailsFile is None:
            self.dbusername = None
            self.dbpassword = None
            self.dbhost = None
            self.dbname = None
        else:
            self.readFile(detailsFile)
            
    def readDetailsLine(self, detailsFile, detailsRe, errorMsg):
        line = detailsFile.readline()
        sys.stderr.write('got line: "%s"\n' % repr(line))
        m = detailsRe.match(line.strip())
        if m is None:
            raise StandardError, errorMsg
        return m.group(1)

    def readFile(self, detailsFile):
        self.dbusername = self.readDetailsLine(detailsFile, self.reDbusername, 'malformed dbusername line')
        self.dbpassword = self.readDetailsLine(detailsFile, self.reDbpassword, 'malformed dbpassword line')
        self.dbhost = self.readDetailsLine(detailsFile, self.reDbhost, 'malformed dbhost line')
        self.dbname = self.readDetailsLine(detailsFile, self.reDbname, 'malformed dbname line')

    def makeConnection(self):
        return mysql.connector.connection.MySQlConnection(user=self.dbusername, password=self.dbpassword, host=self.dbhost, database=self.dbname)

    
