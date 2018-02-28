#!/usr/bin/python

import Tkinter as tk
from Tkinter import *
import Tkinter, Tkconstants, tkFileDialog

from os.path import basename

class BioSankey:

    def __init__(self, master):

        self.master = master
        self.frame = tk.Frame(self.master)
        self.frame.configure()

        n_size=30
        cur_r=1
        
        self.space = Label(self.frame, height=1, text="")
        self.space.grid(row=0,column=0)
        self.help = Text(self.frame, height=4, width=60, wrap=WORD)
        self.help.insert(INSERT, "Host-microbe mapper to analyse dual RNA-seq experiments and determine rRNA genes from NGS data")
        self.help.grid(row=cur_r,column=0, columnspan=3)
        cur_r=cur_r+1

        self.e1A_text = Label(self.frame, width=n_size, text="FASTQ file fwd", anchor='w')
        self.e1A_file = Label(self.frame, width=18, text="", anchor='w', bg = "white", padx=6, borderwidth=2, relief="sunken")
        self.e1A_btn = Button(self.frame, width=5, height=2, bg='lightgrey', text='Select', command=self.degs1A)
        self.e1A_text.grid(row=cur_r,column=0)
        self.e1A_file.grid(row=cur_r,column=1)
        self.e1A_btn.grid(row=cur_r,column=2)
        cur_r=cur_r+1

        self.e1A2_text = Label(self.frame, width=n_size, text="FASTQ file rvs", anchor='w')
        self.e1A2_file = Label(self.frame, width=18, text="", anchor='w', bg = "white", padx=6, borderwidth=2, relief="sunken")
        self.e1A2_btn = Button(self.frame, width=5, height=2, bg='lightgrey', text='Select', command=self.degs1A2)
        self.e1A2_text.grid(row=cur_r,column=0)
        self.e1A2_file.grid(row=cur_r,column=1)
        self.e1A2_btn.grid(row=cur_r,column=2)
        cur_r=cur_r+1
        
        self.e1B_text = Label(self.frame, width=n_size, text="Host ref 1: ", anchor='w')
        self.e1B_file = Label(self.frame, width=18, text="", anchor='w', bg = "white", padx=6, borderwidth=2, relief="sunken")
        self.e1B_btn = Button(self.frame, width=5, height=2, bg='lightgrey', text='Select', command=self.degs1B)
        self.e1B_text.grid(row=cur_r,column=0)
        self.e1B_file.grid(row=cur_r,column=1)
        self.e1B_btn.grid(row=cur_r,column=2)
        cur_r=cur_r+1

        self.e2A_text = Label(self.frame, width=n_size, text="Host ref 2: ", anchor='w')
        self.e2A_file = Label(self.frame, width=18, text="", anchor='w', bg = "white", padx=6, borderwidth=2, relief="sunken")
        self.e2A_btn = Button(self.frame, width=5, height=2, bg='lightgrey', text='Select', command=self.degs2A)
        self.e2A_text.grid(row=cur_r,column=0)
        self.e2A_file.grid(row=cur_r,column=1)
        self.e2A_btn.grid(row=cur_r,column=2)
        cur_r=cur_r+1
        
        self.e2B_text = Label(self.frame, width=n_size, text="Microbe ref 1: ", anchor='w')
        self.e2B_file = Label(self.frame, width=18, text="", anchor='w', bg = "white", padx=6, borderwidth=2, relief="sunken")
        self.e2B_btn = Button(self.frame, width=5, height=2, bg='lightgrey', text='Select', command=self.degs2B)
        self.e2B_text.grid(row=cur_r,column=0)
        self.e2B_file.grid(row=cur_r,column=1)
        self.e2B_btn.grid(row=cur_r,column=2)
        cur_r=cur_r+1


        self.e3A_text = Label(self.frame, width=n_size, text="Host annotation ref1: ", anchor='w')
        self.e3A_file = Label(self.frame, width=18, text="", anchor='w', bg = "white", padx=6, borderwidth=2, relief="sunken")
        self.e3A_btn = Button(self.frame, width=5, height=2, bg='lightgrey', text='Select', command=self.degs3A)
        self.e3A_text.grid(row=cur_r,column=0)
        self.e3A_file.grid(row=cur_r,column=1)
        self.e3A_btn.grid(row=cur_r,column=2)
        cur_r=cur_r+1
        
        self.e3B_text = Label(self.frame, width=n_size, text="Host annotation ref2:", anchor='w')
        self.e3B_file = Label(self.frame, width=18, text="", anchor='w', bg = "white", padx=6, borderwidth=2, relief="sunken")
        self.e3B_btn = Button(self.frame, width=5, height=2, bg='lightgrey', text='Select', command=self.degs3B)
        self.e3B_text.grid(row=cur_r,column=0)
        self.e3B_file.grid(row=cur_r,column=1)
        self.e3B_btn.grid(row=cur_r,column=2)
        cur_r=cur_r+1
        
        self.e2_text = Label(self.frame, width=n_size, text="Microbe annotation ref1", anchor='w')
        self.e2_file = Label(self.frame, width=18, text="", anchor='w', bg = "white", padx=6, borderwidth=2, relief="sunken")
        self.e2_btn = Button(self.frame, width=5, height=2, bg='lightgrey', text='Select', command=self.degsE2)
        self.e2_text.grid(row=cur_r,column=0)
        self.e2_file.grid(row=cur_r,column=1)
        self.e2_btn.grid(row=cur_r,column=2)
        cur_r=cur_r+1

        self.e3_text = Label(self.frame, width=n_size, text="Testrun: ", anchor='w')
        self.e3_btn = Scale(self.frame, from_=0, to=100, orient=HORIZONTAL)
        self.e3_text.grid(row=cur_r,column=0)
        self.e3_btn.grid(row=cur_r,column=1)
        cur_r=cur_r+1
        
        self.CheckVar1=IntVar()
        self.e4_text = Label(self.frame, width=n_size, text="Assembly of 16S rRNA genes: ", anchor='w')
        self.e4_btn = Checkbutton(self.frame, text = "", variable = self.CheckVar1, onvalue = 1, offvalue = 0, height=5, width = 20)
        self.e4_text.grid(row=cur_r,column=0)
        self.e4_btn.grid(row=cur_r,column=1)
        cur_r=cur_r+1
                 
        self.startBiosankey = tk.Button(self.frame, height=2, text = 'Start', width = 5, command = self.start_biosankey, bg="blue")
        self.startBiosankey.grid(row=cur_r,column=1)
        cur_r=cur_r+1

        self.f_degs=None
        self.frame.pack()

    # -------------------------------------------------------------------------------

    def degs1A(self):
        self.f_degs = tkFileDialog.askopenfilename()
        if(self.f_degs): 
            self.e1A_btn.config(text='Reset', command=self.reset_degs1A)
            self.e1A_file.config(text=basename(self.f_degs))

    def degs1A2(self):
        self.f_degs = tkFileDialog.askopenfilename()
        if(self.f_degs):
            self.e1A2_btn.config(text='Reset', command=self.reset_degs1A2)
            self.e1A2_file.config(text=basename(self.f_degs))

    def degs1B(self):
        self.f_degs = tkFileDialog.askopenfilename()
        if(self.f_degs):
            self.e1B_btn.config(text='Reset', command=self.reset_degs1B)
            self.e1B_file.config(text=basename(self.f_degs))

    def degs2A(self):
        self.f_degs = tkFileDialog.askopenfilename()
        if(self.f_degs):
            self.e2A_btn.config(text='Reset', command=self.reset_degs2A)
            self.e2A_file.config(text=basename(self.f_degs))

    def degs2B(self):
        self.f_degs = tkFileDialog.askopenfilename()
        if(self.f_degs):
            self.e2B_btn.config(text='Reset', command=self.reset_degs2B)
            self.e2B_file.config(text=basename(self.f_degs))

    def degs3A(self):
        self.f_degs = tkFileDialog.askopenfilename()
        if(self.f_degs):
            self.e3A_btn.config(text='Reset', command=self.reset_degs3A)
            self.e3A_file.config(text=basename(self.f_degs))

    def degs3B(self):
        self.f_degs = tkFileDialog.askopenfilename()
        if(self.f_degs):
            self.e3B_btn.config(text='Reset', command=self.reset_degs3B)
            self.e3B_file.config(text=basename(self.f_degs))

    def degsE2(self):
        self.f_degs = tkFileDialog.askopenfilename()
        if(self.f_degs):

            self.e2_btn.config(text='Reset', command=self.reset_degsE2)
            self.e2_file.config(text=basename(self.f_degs))

    # -------------------------------------------------------------------------------

    def reset_degs1A(self):
        self.e1A_file.config(text="")
        self.e1A_btn.configure(text='Select', command=self.degs1A)

    def reset_degs1A2(self):
        self.e1A2_file.config(text="")
        self.e1A2_btn.configure(text='Select', command=self.degs1A2)

    def reset_degs1B(self):
        self.e1B_file.config(text="")
        self.e1B_btn.configure(text='Select', command=self.degs1B)

    def reset_degs2A(self):
        self.e2A_file.config(text="")
        self.e2A_btn.configure(text='Select', command=self.degs2A)

    def reset_degs2B(self):
        self.e2B_file.config(text="")
        self.e2B_btn.configure(text='Select', command=self.degs2B)

    def reset_degs3A(self):
        self.e3A_file.config(text="")
        self.e3A_btn.configure(text='Select', command=self.degs3A)

    def reset_degs3B(self):
        self.e3B_file.config(text="")
        self.e3B_btn.configure(text='Select', command=self.degs3B)

    def reset_degsE2(self):
        self.e2_file.config(text="")
        self.e2_btn.configure(text='Select', command=self.degsE2)

    # -------------------------------------------------------------------------------

    def close_windows(self):
        self.master.destroy()

    def start_biosankey(self):
	e1A=self.e1A_file["text"]
	e1A2=self.e1A2_file["text"]
        e1B=self.e1B_file["text"]
        e2A=self.e2A_file["text"]
        e2B=self.e2B_file["text"]
        e3A=self.e3A_file["text"]
        e3B=self.e3B_file["text"]
	e3=self.e3_btn.get()
	e4=self.CheckVar1.get()

	if(len(e1A)>0 and len(e1A2)>0 and len(e1B)>0 and len(e2B)>0 and len(e3A)>0):
		print "FASTQ FWD FASTQ: "+e1A
		print "FASTQ RVS FASTQ: "+e1A2
		print "Host ref 1: "+e1B
		print "Host ref 2: "+e2A
		print "Microbe ref 1: "+e2B
		print "Host anno 1: "+e3A
		print "Host anno 2: "+e3B
		print "Testrun: "+str(e3)
		print "16S Assembly: "+str(e4)
		print "-----------------------------------------------------"

		template="""#!/bin/bash
module load nano
module load java

PATH=~/.local/bin/:$PATH	
PATH=bin/Trimmomatic-0.36:$PATH
PATH=bin/sortmerna-2.1b/:$PATH
PATH=/naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/bin/star/2.5.3a/bin/Linux_x86_64/:$PATH
PATH=/naslx/projects_mpiio/pr74ma/ge34juq2/HOST_MICROBE_MAPPER/host_pathogen_mapping-master/bin/tophat/2.1.1:$PATH
PATH=bin/bwa-0.7.17/:$PATH
PATH=bin/FastQC/:$PATH
PATH=bin/sortmerna-2.1b/scripts/:$PATH

./data/host_pathogen_mapping.sh -F __e1A__ -R __e1A2__ -P __e2B__ -C 1 -X output -H __e1B__ -I __e2A__ -J __e3A__ -K __e3B__ -V __e3__
		"""
	
		template=template.replace("__e1A__", e1A)
                template=template.replace("__e1A2__",e1A2)
                template=template.replace("__e1B__", e1B)
                template=template.replace("__e2A__", e2A)
                template=template.replace("__e2B__", e2B)
                template=template.replace("__e3A__", e3A)
                template=template.replace("__e3B__", e3B)
                template=template.replace("__e3__",  str(e3))
                template=template.replace("__e4__",  str(e4))

		fw=file("pipeline_v2.sh","w")
		fw.write(template)
		fw.close()

	else:
		print "INFO\tValues are missing"
        
def main(): 
    root = tk.Tk()
    root.title("HostMicrobeMapper - version 0.1")
    app = BioSankey(root)
    root.configure(background='white')
    root.geometry('{}x{}'.format(500, 550))
    root.mainloop()

if __name__ == '__main__':
    main()

    
    


