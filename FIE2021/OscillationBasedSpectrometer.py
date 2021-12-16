import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import WolframSolver as WS
import scipy.io as sio
import symlib as sl

class OBS:
    def __init__(self, BioZ, CalZ, Comp, G, Q, f_list):
        self.BZ = BioZ
        self.CalZ = CalZ
        self.C = Comp
        self.G = G
        self.Q = Q
        self.f_list = f_list
        self.i = 0
        self.WF = WS.Wolfram()

    def runJob(self):
        self.i = 0
        if (self.CalZ is not None):
            jobdict = {
                'nFreq': len(self.f_list),
                'cal': True,
                'joblist': []
            }
            while(self.i < len(self.f_list)):
                BPF_ = sl.BPF(Q=self.Q, w0=self.f_list[self.i], G0=self.G)
                OBT_ = sl.OBT(BPF_, self.CalZ, self.C)
                eq, sym = OBT_.getEqusym()
                job = {
                    "fr": self.f_list[self.i],
                    "equ": eq,
                    "sym": sym,
                }
                jobdict['joblist'].append(job)
                self.i += 1
            self.WF.createJob(jobdict)
            self.WF.execute()

        self.i = 0
        jobdict = {
            'nFreq': len(self.f_list),
            'cal': False,
            'joblist': []
        }
        while (self.i < len(self.f_list)):
            BPF_ = sl.BPF(Q=self.Q, w0=self.f_list[self.i], G0=self.G)
            OBT_ = sl.OBT(BPF_, self.BZ, self.C)
            eq, sym = OBT_.getEqusym()
            job = {
                "fr": self.f_list[self.i],
                "equ": eq,
                "sym": sym,
            }
            jobdict['joblist'].append(job)
            self.i += 1

        self.WF.createJob(jobdict)
        self.WF.execute()


    def saveMat(self, filename, key):
        wolfram_dir = (os.path.abspath(os.getcwd())+"\Wolfram\\").replace("\\", "\\\\")
        output_dir = (os.path.abspath(os.getcwd())+"\matfiles\\").replace("\\", "\\\\")
        output_file = open(output_dir + filename, "wb")

        F = []
        AC = []
        A = []
        WC = []
        W = []

        try:
            fileC = open(wolfram_dir+"ResultCAL.csv")
            cal_csv = csv.reader(fileC)
            for row in cal_csv:
                w = float(row[1])
                a = float(row[2])
                WC.append(w)
                AC.append(a)
            fileC.close()
            os.remove(wolfram_dir+"ResultCAL.csv")
        except Exception  as e:
            pass

        try:
            file = open(wolfram_dir + "Result.csv")
            res_csv = csv.reader(file)
            for row in res_csv:
                fr = float(row[0])
                try:
                    w = float(row[1])
                    a = float(row[2])
                except:
                    w = -1
                    a = -1
                F.append(fr)
                W.append(w)
                A.append(a)
            file.close()
            os.remove(wolfram_dir + "Result.csv")
        except Exception  as e:
            pass

        mdict = {"F": F, "A": A, "AC": AC, "W": W, "WC": WC, "Cal": len(AC) >0}
        data = {key: mdict}
        sio.savemat(output_file, data)


    def appendMat(self, filename, key):
        wolfram_dir = (os.path.abspath(os.getcwd())+"\Wolfram\\").replace("\\", "\\\\")
        output_dir = (os.path.abspath(os.getcwd())+"\matfiles\\").replace("\\", "\\\\")

        if(not os.path.isfile(output_dir + filename)):
            output_file = open(output_dir + filename, "wb")
            data = {}
        else:
            data = sio.loadmat(output_dir + filename)
            output_file = open(output_dir + filename, "wb")

        F = []
        AC = []
        A = []
        WC = []
        W = []

        try:
            fileC = open(wolfram_dir+"ResultCAL.csv")
            cal_csv = csv.reader(fileC)
            for row in cal_csv:
                w = float(row[1])
                a = float(row[2])
                WC.append(w)
                AC.append(a)
            fileC.close()
            os.remove(wolfram_dir+"ResultCAL.csv")
        except Exception  as e:
            pass

        try:
            file = open(wolfram_dir + "Result.csv")
            res_csv = csv.reader(file)
            for row in res_csv:
                fr = float(row[0])
                try:
                    w = float(row[1])
                    a = float(row[2])
                except:
                    w = -1
                    a = -1
                F.append(fr)
                W.append(w)
                A.append(a)
            file.close()
            os.remove(wolfram_dir + "Result.csv")
        except Exception  as e:
            pass

        mdict = {"F": F, "A": A, "AC": AC, "W": W, "WC": WC, "Cal": len(AC) >0}
        data[key] = mdict
        sio.savemat(output_file, data)

