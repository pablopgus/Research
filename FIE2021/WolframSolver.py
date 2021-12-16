import os
import re
import contextlib
import sys, os

class Wolfram:
    def __init__(self):
        self.dir_path = (os.path.abspath(os.getcwd())+"\Wolfram\\").replace("\\", "\\\\")
        self.p = re.compile('[A,W] -> [0-9.]+')
        self.d = re.compile('[0-9.]+')

    def execute(self):
        with contextlib.redirect_stdout(None):
            os.system("math -script "+self.dir_path+"SolveEQ.wls")

    def createJob(self, job_dict):
        f = open(self.dir_path+"SolveEQ.wls", "w")
        header = "# !/usr/bin/env wolframscript\n(* ::Package:: *)\n\n"
        if (job_dict['cal']):
            header += "file = \"" + self.dir_path + "ResultCAL.csv\"\n\n"
        else:
            header += "file = \"" + self.dir_path + "Result.csv\"\n\n"

        #header += "DeleteFile[file]\n\n"

        header += "nFreq = "+str(job_dict["nFreq"])+";\n DATA = Array[0.0, {nFreq, 3}]\n\n"

        body = ""
        i = 1;
        for job in job_dict["joblist"]:
            t_str = self.createS_str(job['equ'], job['sym'])
            body += t_str
            body += "RES = S[[All, All, 2]]\nw = RES[[1]][[1]]\na = RES[[1]][[2]]\n\n"
            body += "DATA[["+str(i)+",1]] = "+str(job['fr'])+"\nDATA[["+str(i)+",2]] = w\nDATA[["+str(i)+",3]] = a\n\n"
            i += 1

        footer = "Export[file,DATA]\n"
        f.write(header + body + footer)


    def createS_str(self, equations, symbols):
        eq = "S = Solve[{ \n"
        for e in equations:
            eq += "\t\t" + str(e).replace("**", "^").replace("e", "*10^") + " == 0, \n"
        eq += "\t\tW > 0, A > 0}, \n\t\t{"

        i = 0
        while i < len(symbols):
            if (i < (len(symbols) - 1)):
                eq += str(symbols[i]) + ","
            else:
                eq += str(symbols[i])
            i += 1
        eq += "}, Reals]\n\n"
        return eq

    def createScript(self, equations, symbols):
        f = open(self.dir_path+"SolveEQ.wls", "w")
        header = "# !/usr/bin/env wolframscript\n(* ::Package:: *)\n\n"
        header += "file = \""+self.dir_path+"Result.res\"\n\n"
        header += "DeleteFile[file]\n\n"

        # A^2 + W = 0, W + 2 = 0}, {A, W}

        eq = "S = Solve[{ \n"
        for e in equations:
            eq += "\t\t"+str(e).replace("**", "^").replace("e", "*10^") + " == 0, \n"
        eq += "\t\tW > 0, A > 0}, \n\t\t{"

        i = 0
        while i < len(symbols):
            if(i < (len(symbols)-1) ):
                eq += str(symbols[i]) + ","
            else:
                eq += str(symbols[i])
            i += 1
        eq += "}, Reals]\n\n"

        footer = "Save[file, S]\n"
        f.write(header + eq + footer)

    def getResults(self):
        f = open(self.dir_path+"Result.res", "r")
        file_res = f.read().replace("\n", "")
        matchs = self.p.findall(file_res)
        try:
            for m in matchs:
                v = float(self.d.findall(m)[0])
                if(m[0] == 'A'):
                    A = v
                elif(m[0] == 'W'):
                    F = v
            return A,F
        except:
            return 0,0
