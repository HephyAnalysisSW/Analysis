"""
Script to create Hephy/CERN tokens needed to run on condor from lxplus
"""

import glob
import subprocess as sp
import shutil
import shlex
import json
import os

tokenPath = os.environ["HOME"] +  "/private/kerberos"
def main():
    prepareTokens()
    useToken("hephy")

def checkTokens():
        prepareTokens()
        return True


def prepareTokens():
        on_system = getSystem()

        p = sp.Popen( shlex.split("klist"), stdout=sp.PIPE, stderr=sp.PIPE )
        (out,err) =  p.communicate()    
        for line in out.splitlines():

            if "Ticket cache" in line:
                token = line.split("FILE:")[1]

        if not os.path.exists(tokenPath):
            os.mkdir(tokenPath)
        if not os.path.exists(tokenPath+"/krb5_token_"+on_system):
            shutil.copyfile(token, tokenPath+"/krb5_token_"+on_system )

        if on_system == "cern.ch":
            if os.path.exists(tokenPath+"/username"):
                with open(tokenPath+"/username","r") as FSO:
                    user = FSO.read()
            else:
                user = input("Give username for cell hephy.at: ")
                save = input("You want to save '{0}' as your username for later? (Y/n)".format(user))
                if save != "n":
                    with open(tokenPath+"/username","w") as FSO:
                        FSO.write(user)

            p = sp.Popen( shlex.split("kinit {0}@HEPHY.AT -c %s/krb5_token_hephy.at".format(user)%tokenPath) )
            p.communicate()
            os.environ["KRB5CCNAME"] = tokenPath + "/krb5_token_hephy.at"
            os.system("aklog -d hephy.at")
            os.environ["KRB5CCNAME"] = tokenPath + "/krb5_token_cern.ch"


def useToken(cell):

    if cell == "hephy":
        os.environ["KRB5CCNAME"] = tokenPath + "/krb5_token_hephy.at"
    if cell == "cern":
        os.environ["KRB5CCNAME"] = tokenPath + "/krb5_token_cern.ch"



def checkProxy():

        p = sp.Popen( shlex.split("voms-proxy-info"), stdout=sp.PIPE, stderr=sp.PIPE )
        (out,err) =  p.communicate()

        if err:
            print(err)
            return False

        info = {"time":-1,"path":""}
        for line in out.splitlines():
            if "timeleft" in line: info['time'] = int( line.replace(" ","").split("timeleft:")[1].replace(":","") )
            if "path" in line: info["path"] = line.split(":")[1].replace(" ","")

        if info['time'] > 0:
            if not os.path.exists("proxy"):
                os.mkdir("proxy")
            shutil.copyfile(info["path"], "proxy/x509_proxy")
            return True
        else:
            print("Proxy expired! Get a new one with 'voms-proxy-init --voms cms'")
            return False

def getHeplxPublicFolder():

    user = os.environ["USER"]
    first = user[0]
    second = user[:8]

    return glob.glob( "/afs/hephy.at//user/{0}/{1}*/public".format(first,second) )[0]

def getSystem(inverse = False):
    host = os.environ["HOSTNAME"]
    if "heplx" in host: 
      if inverse: return "cern.ch"
      else: return "hephy.at"

    elif "lxplus" in host: 
      if inverse: return "hephy.at"
      else: return "cern.ch"

if __name__ == '__main__':
    main()
