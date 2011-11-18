# ---------------------------------------------------------------
#
# Python tools for GTMK
# Simon King, 2004
# beta version, not for redistribution
#
# ---------------------------------------------------------------

import os
import time
import re
import threading
import popen2
import signal
import pwd

import Tkinter
import Pmw

# ---------------------------------------------------------------

def prepend(l,i):
    # dumb and probably slow
    l.reverse()
    l.append(i)
    l.reverse()


class myPmwCounter(Pmw.Counter):

    myCallback=None
    def __init__(self, parent = None, callback=None, **kw):
        self.myCallback=callback
	apply(Pmw.Counter.__init__, (self, parent), kw)
        #print "   made myPmwCounter",self.myCallback,self
        
    # just override the function called after an arrow button is released
    def _stopCounting(self, event = None):
        Pmw.Counter._stopCounting(self,event)
        #print "stop; callback=",self.myCallback
        #print "self=",self
        if self.myCallback:
            apply(self.myCallback)
        
        
class Displayer(threading.Thread):

    # multi-purpose display class
    # consists of tabs, each for a specific purpose
    # can add new tabs for new tasks, like viewing DTs, editing DBNs (!?)

    import MultiListbox as MLB

    master=None
    agenda=None
    list=None
    button=update_button=None
    q=None
    interval=int(30) # seconds
    guiMainThread=None
    gui_initialised=False
    text=[]
    agenda_update_buttons=[]

    def __init__(self,
                 agendas=None,
                 lock=None, tablist=[],
                 target=None, name='gmtkPT_Displayer',
                 parameters={}):

        self.parameters=parameters
        self.thread_lock=lock
        if self.thread_lock:
            self.thread_lock.acquire()
        else:
            print "Warning: GUI starting without thread locking"

        self.re_agenda=re.compile(r'^agenda\:(\w+)')
            
        print "agendas=",agendas

        self.agendas=agendas # the agenda threads
        self.tablist=tablist

        if self.agendas:
            for name in self.agendas.keys():
                self.tablist.append(name)
        
        self.gui_status="notReady"

        self.actions_page=None
        self.update_page=None
        self.monitor_page=None
        self.train_page=None
        self.makeInitialParameters_page=None

        # those Tk widgets we need to retain so we can update them
        # stored as dictionaries keyed by tab name
        self.agenda_pages={}
        self.progress_bars={}
        self.remaining_time={}
        self.output_group_label_texts={}
        self.lists={}
        self.bestboxes={}
        self.cmdboxes={}
        self.textboxes={}
            
        # get some info on the current process
        self.local_process_monitor_users=None
        self.this_pid=os.getpid()
        self.this_uid=os.getuid()
        try:
            self.this_username=os.getlogin()
        except OSError:
            self.this_username=pwd.getpwuid(os.getuid())[0]

        if not target:
            target=self.start_mainloop

        # start this thread
        threading.Thread.__init__(self, target=target, name=name, args=())
        self.setDaemon(False)
        self.start()

    def __delete__(self):
        print "GUI object deleted"
        self.die()


    def set_interval(self,i):

        # update interval controlled by Tkinter timer
        if self.timer_id:
            self.update_button.after_cancel(self.timer_id)

        self.interval=int(i)
        self.timer_id=self.update_button.after(self.interval*1000, self.start_timer)

    def get_interval(self):
        return self.interval

    def start_mainloop(self):
        #try:
            self.initialise_display()
            self.gui_initialised=True
            self.gui_status="initialising"


            print "starting main loop, thread is",threading.currentThread()
            self.guiMainThread=threading.currentThread()

            self.start_timer()
            self.gui_status="running"
            
            if self.thread_lock:
                self.thread_lock.release()

            self.master.mainloop()

            self.gui_status="finished"
            self.gui_initialised=False
            #self.die()

        #except Tkinter.TclError:
        #    print "Got a Tcl error during intitialisation"
        #    self.gui_initialised=False
        #    self.gui_status="initialisationFailed"
            #self.die()


    def initialise_display(self,title='gmtkPT Display'):
        self.master = Tkinter.Tk() # main window
        self.master.config()
        self.master.title(title)

        # need this variable even if there is no monitor tab
        self.local_process_monitor_users = Tkinter.StringVar()


        # --------------------------------------------------------------
	# Create a Pmw NoteBook (a tabbed set of pages)
        self.notebook = Pmw.NoteBook(self.master)
        self.notebook.pack(fill = 'both', expand = 1, padx = 0, pady = 0)

        for t in self.tablist:
            self.create_tab(t)


    def create_tab(self,tabname):

        if tabname =="actions":
            self.create_actions_tab()
        elif tabname =="settings":
            self.create_settings_tab()
        elif tabname =="updateSettings":
            self.create_updateSettings_tab()
        elif tabname =="monitor":
            self.create_monitor_tab()
        elif self.re_agenda.match(tabname):
            self.create_agenda_tab(name=tabname)
        elif tabname =="train":
            self.create_train_tab()
        elif tabname =="makeInitialParameters":
            self.create_makeInitialParameters_tab()
        else:
            raise ValueError, "Tried to create unknown tab: "+tabname



    def create_actions_tab(self):
        self.actions_page = self.notebook.add('Actions')
        self.notebook.tab('Actions').focus_set()

        self.actions_group = Pmw.Group(self.actions_page, tag_text = 'Actions')
        self.actions_group.pack(fill = 'both', expand = 0, padx = 0, pady = 0)

        self.button = Tkinter.Button(self.actions_group.interior(),
                                          text="QUIT", fg="red", command=self.quit_button_pressed)
        self.button.grid(column=0,row=0)

        self.update_button = Tkinter.Button(self.actions_group.interior(),
                                                 text="Update", command=self.update)
        # not packed - so actually invisible
        #self.update_button.grid(column=1,row=0)

    def create_updateSettings_tab(self):
        self.update_page = self.notebook.add('Update settings')
        self.notebook.tab('Update settings').focus_set()

        self.interval_slider = Tkinter.Scale(self.update_page, label="Update interval (s)",
                                                  from_=1, to=600, orient=Tkinter.HORIZONTAL,
                                                  length=200, command=self.set_interval)
        self.interval_slider.set(self.interval)
        self.interval_slider.grid(column=2, row=0)


    def create_agenda_tab(self,name="Agenda"):
        self.agenda_pages[name] = self.notebook.add(name)
        self.notebook.tab(name).focus_set()

        self.agenda_update_button = Tkinter.Button(self.agenda_pages[name],
                                                 text="Update", command=self.update)
        self.agenda_update_button.pack(side=Tkinter.TOP)
        self.agenda_update_buttons.append(self.agenda_update_button)

        self.progress_bars[name] = ProgressBar(self.agenda_pages[name],
                                               value=0)
        
        self.progress_bars[name].frame.pack(side=Tkinter.TOP)

        self.remaining_time[name] = Tkinter.StringVar()
        self.remaining_time[name].set('Remaining time unknown')
        this_text = Tkinter.Label(self.agenda_pages[name],
                                  textvariable=self.remaining_time[name])
        this_text.pack(side=Tkinter.TOP)

        agenda_page_panes = Pmw.PanedWidget(self.agenda_pages[name],
                                            orient='vertical',
                                            hull_borderwidth = 1,
                                            hull_relief = 'sunken')
        
        agenda_page_panes.pack(expand = 1, fill='both')
        
        pane = agenda_page_panes.add("bestPane", min = .1, size = .2)
        best_group = Pmw.Group(pane, tag_text = 'Best so far')
        best_group.pack(fill = 'both', expand = 1, padx = 0, pady = 0)
        
        self.bestboxes[name] = Tkinter.Text(best_group.interior(), height=10,width=200)
        bestbox_scrollbar = Tkinter.Scrollbar(best_group.interior())
        bestbox_scrollbar.config(command=self.bestboxes[name].yview)
        self.bestboxes[name].config(yscrollcommand=bestbox_scrollbar.set)
        bestbox_scrollbar.pack(side=Tkinter.RIGHT,fill=Tkinter.Y)
        self.bestboxes[name].pack(side=Tkinter.LEFT,fill=Tkinter.X)
        
        
        pane = agenda_page_panes.add("outputPane", min = .1, size = .4)
        output_group = Pmw.Group(pane, tag_text='Job')
        output_group.pack(fill = 'both', expand = 1, padx = 0, pady = 0)
        
        self.output_group_label_texts[name] = Tkinter.StringVar()
        self.output_group_label_texts[name].set('Double click a job to see output here')
        output_group_label = Tkinter.Label(output_group.interior(),
                                           textvariable=self.output_group_label_texts[name])
        
        output_group_label.pack()
        
        output_group_cmd = Pmw.Group(output_group.interior(), tag_text='Command')
        output_group_cmd.pack(fill = 'both', expand = 1, padx = 0, pady = 0)
        
        self.cmdboxes[name] = Tkinter.Text(output_group_cmd.interior())
        cmdbox_scrollbar = Tkinter.Scrollbar(output_group_cmd.interior())
        cmdbox_scrollbar.config(command=self.cmdboxes[name].yview)
        self.cmdboxes[name].config(yscrollcommand=cmdbox_scrollbar.set,height=10,width=200)
        cmdbox_scrollbar.pack(side=Tkinter.RIGHT,fill=Tkinter.Y)
        self.cmdboxes[name].pack(side=Tkinter.LEFT,fill=Tkinter.X)
        
        output_group_text = Pmw.Group(output_group.interior(), tag_text='Output')
        output_group_text.pack(fill = 'both', expand = 1, padx = 0, pady = 0)
        
        self.textboxes[name] = Tkinter.Text(output_group_text.interior())
        textbox_scrollbar = Tkinter.Scrollbar(output_group_text.interior())
        textbox_scrollbar.config(command=self.textboxes[name].yview)
        self.textboxes[name].config(yscrollcommand=textbox_scrollbar.set,height=10,width=200)
        textbox_scrollbar.pack(side=Tkinter.RIGHT,fill=Tkinter.Y)
        self.textboxes[name].pack(side=Tkinter.LEFT,fill=Tkinter.X)
            
        pane = agenda_page_panes.add("jobsPane", min = .1, size = .4)
        jobs_group = Pmw.Group(pane, tag_text = 'Jobs')
        jobs_group.pack(fill = 'both', expand = 1, padx = 0, pady = 0)
        
        self.lists[name] = self.MLB.MultiListbox(jobs_group.interior(),
                                                 ( ('Job', 7), ('S', 10), ('Dep', 6),
                                                   ('Part', 12), ('Node', 10),
                                                   ('Command', 16), ('Output', 55)  ) )
        self.lists[name].pack()
        
        self.lists[name].bind("<Double-Button-1>", self.select_one)



    def create_monitor_tab(self):

        self.monitor_page = self.notebook.add('Monitor')
        self.notebook.tab('Monitor').focus_set()

        self.monitor_update_button = Tkinter.Button(self.monitor_page,
                                                 text="Update", command=self.update)
        self.monitor_update_button.pack()

        
        self.local_group = Pmw.Group(self.monitor_page, tag_text = 'Local processes')
        self.local_group.pack(fill = 'both', expand = 1, padx = 10, pady = 10)
        
        self.local_b1=Tkinter.Radiobutton(self.local_group.interior(),
                                               text="Children of this process ("+self.this_pid.__repr__()+")",
                                               variable=self.local_process_monitor_users, value="children")
        b2=Tkinter.Radiobutton(self.local_group.interior(), text="All of "+self.this_username+"'s processes",
                                    variable=self.local_process_monitor_users, value="self")
        b3=Tkinter.Radiobutton(self.local_group.interior(), text="All processes",
                                    variable=self.local_process_monitor_users, value="all")
        
        self.local_b1.select()
        self.local_b1.pack(anchor=Tkinter.W)
        b2.pack(anchor=Tkinter.W)
        b3.pack(anchor=Tkinter.W)
        
        self.kill_local_button = Tkinter.Button(self.local_group.interior(),
                                                     text="Kill the children!", command=self.kill_local)
        
        self.kill_local_button.pack()
        
        self.localbox = Tkinter.Text(self.local_group.interior(),width=200)
        self.localbox_scrollbar = Tkinter.Scrollbar(self.local_group.interior())
        self.localbox_scrollbar.config(command=self.localbox.yview)
        self.localbox.config(yscrollcommand=self.localbox_scrollbar.set)
        self.localbox_scrollbar.pack(side=Tkinter.RIGHT,fill=Tkinter.Y)
        self.localbox.pack(side=Tkinter.LEFT,fill=Tkinter.X)
        
        self.remote_group = Pmw.Group(self.monitor_page, tag_text = 'Remote processes')
        self.remote_group.pack(fill = 'both', expand = 1, padx = 10, pady = 10)
        
        self.remote_process_monitor_users = Tkinter.StringVar()
        b2=Tkinter.Radiobutton(self.remote_group.interior(), text="All of "+self.this_username+"'s processes",
                                    variable=self.remote_process_monitor_users, value="self")
        b3=Tkinter.Radiobutton(self.remote_group.interior(), text="All processes",
                                    variable=self.remote_process_monitor_users, value="all")
        
        b2.select()
        b2.pack(anchor=Tkinter.W)
        b3.pack(anchor=Tkinter.W)
        
        #self.kill_remote_button = Tkinter.Button(self.remote_group.interior(),
        #                                             text="Kill all remote jobs!", command=self.kill_remote)
        
        #self.kill_remote_button.pack()
        
        self.remote_group_label_text = Tkinter.StringVar()
        self.remote_group_label_text.set('Using cctrl')
        self.remote_group_label = Tkinter.Label(self.remote_group.interior(),
                                                     textvariable=self.remote_group_label_text)
        self.remote_group_label.pack()
            
        self.remotebox = Tkinter.Text(self.remote_group.interior(),width=200)
        self.remotebox_scrollbar = Tkinter.Scrollbar(self.remote_group.interior())
        self.remotebox_scrollbar.config(command=self.remotebox.yview)
        self.remotebox.config(yscrollcommand=self.remotebox_scrollbar.set)
        self.remotebox_scrollbar.pack(side=Tkinter.RIGHT,fill=Tkinter.Y)
        self.remotebox.pack(side=Tkinter.LEFT,fill=Tkinter.X)
        

    def create_train_tab(self):
        self.train_page = self.notebook.add('Training')
        self.notebook.tab('Training').focus_set()
        # not yet implemented


    def custom_validator_function(self, min, max, text):
        print min,max,text
        i=int(text)
        print "validating DCPT"
        return ((i >= 0) and (i < 100))


    def create_makeInitialParameters_tab(self):
        self.makeInitialParameters_page = self.notebook.add('Initial parameter creator')
        self.notebook.tab('Initial parameter creator').focus_set()

        # parameters that can be made are currently:
        # Dense CPTs (for transition matrices)
        # Gaussian mixture distributions (Gaussians, DPMFs for weights)

        self.makeInitialParameters_page_panes = Pmw.PanedWidget(self.makeInitialParameters_page,
                                                                     orient='vertical',
                                                                     hull_borderwidth = 1,
                                                                     hull_relief = 'sunken')
        
        self.makeInitialParameters_page_panes.pack(expand = 1, fill='both')

        pane = self.makeInitialParameters_page_panes.add("pane1", min = .1, size = .1)
        self.makeInitialParameters_actions_group = Pmw.Group(pane, tag_text = 'Actions')
        self.makeInitialParameters_actions_group.pack(fill = 'both', expand = 1, padx = 10, pady = 10)


        self.makeInitialParameters_save_button = Tkinter.Button(self.makeInitialParameters_actions_group.interior(),
                                                                text="Save",
                                                                command=self.save)
        self.makeInitialParameters_save_button.pack()


        pane = self.makeInitialParameters_page_panes.add("pane2", min = .1, size = .5)
        self.makeInitialParameters_preview_group = Pmw.Group(pane, tag_text = 'Preview')
        self.makeInitialParameters_preview_group.pack(fill = 'both', expand = 1, padx = 0, pady = 0)

        self.preview_box = Tkinter.Text(self.makeInitialParameters_preview_group.interior())
        self.preview_box_scrollbar = Tkinter.Scrollbar(self.makeInitialParameters_preview_group.interior())
        self.preview_box_scrollbar.config(command=self.preview_box.yview)
        self.preview_box.config(yscrollcommand=self.preview_box_scrollbar.set,height=30,width=200,background='#FFFFFF')
        self.preview_box_scrollbar.pack(side=Tkinter.RIGHT,fill=Tkinter.Y)
        self.preview_box.pack(side=Tkinter.LEFT,fill=Tkinter.BOTH)

        pane = self.makeInitialParameters_page_panes.add("pane3", min = .1, size = .2)
        self.DCPT_group = Pmw.ScrolledFrame(pane) # , tag_text = 'Dense CPTs')
        self.DCPT_group.pack(fill = 'both', expand = 1, padx = 0, pady = 0)

        self.num_dcpts = myPmwCounter(self.DCPT_group.interior(),
                                      callback=self.make_dcpts,
                                      datatype='integer',
                                      labelpos = 'w',
                                      label_text = 'How many DCPTs?',
                                      orient = 'horizontal',
                                      entry_width = 2,
                                      entryfield_value = 0,
                                      entryfield_command=self.make_dcpts,
                                      entryfield_validate = {'validator' : 'integer',
                                                             'min' : 0, 'max' : 99})

        
        #self.dcpt_button = Tkinter.Button(self.DCPT_group.interior(),
        #                                       text="Update",
        #                                       command=self.make_dcpts)
        
        self.num_dcpts.pack()
        #self.dcpt_button.pack()
        self.dcpts=[]


        pane = self.makeInitialParameters_page_panes.add("pane4", min = .1, size = .2)

        self.gaussian_group = Pmw.ScrolledFrame(pane)
        self.gaussian_group.pack(fill = 'both', expand = 1, padx = 0, pady = 0)

        #self.gaussian_group = Pmw.Group(self.gaussian_group_container,
        #                                tag_text = 'Gaussian collections')
        #self.gaussian_group.pack(fill = 'both', expand = 1, padx = 10, pady = 10)


        #self.gaussian_scrollbar = Tkinter.Scrollbar(self.gaussian_group.interior())
        #self.gaussian_scrollbar.config(command=self.gaussian_group.yview)
        #self.xxxx.config(yscrollcommand=self.gaussian_scrollbar.set)
        #self.gaussian_scrollbar.pack(side=Tkinter.RIGHT,fill=Tkinter.Y)

        self.num_gaussian_collections = myPmwCounter(self.gaussian_group.interior(),
                                                     callback=self.make_gaussians,
                                                     datatype='integer',
                                                     labelpos = 'w',
                                                     label_text = 'How many collections of Gaussian mixture distributions?',
                                                     orient = 'horizontal',
                                                     entry_width = 2,
                                                     entryfield_value = 0,
                                                     entryfield_command=self.make_gaussians,
                                                     entryfield_validate = {'validator' : 'integer',
                                                                            'min' : 0, 'max' : 99})
        #self.gaussian_button = Tkinter.Button(self.gaussian_group.interior(),
        #                                      text="Update",
        #                                      command=self.make_gaussians)
        
        self.num_gaussian_collections.pack()
        #self.gaussian_button.pack()
        self.gaussian_collections=[]


    def make_dcpts(self):
        print "making DCPTs"
        if not self.DCPT_group:
            return

        while int(self.num_dcpts.get()) < len(self.dcpts):
            print "forgetting down to",int(self.num_dcpts.get())
            self.dcpts[len(self.dcpts)-1].pack_forget()
            del self.dcpts[len(self.dcpts)-1]

        while int(self.num_dcpts.get()) > len(self.dcpts):
            print "adding for final length of",int(self.num_dcpts.get()) 
            self.dcpts.append(DcptField(self.DCPT_group.interior(),
                                        tag_text = 'DCPT '+len(self.dcpts).__repr__()))
            self.dcpts[len(self.dcpts)-1].create(name='dcpt'+len(self.dcpts).__repr__(),
                                                 callback=self.makeInitialParameters_preview)
            self.dcpts[len(self.dcpts)-1].pack(fill=Tkinter.X)

        self.makeInitialParameters_preview()
        
    def make_gaussians(self):
        print "making gaussian mixtures"
        if not self.gaussian_group:
            return

        while int(self.num_gaussian_collections.get()) < len(self.gaussian_collections):
            print "forgetting down to",int(self.num_gaussian_collections.get())
            self.gaussian_collections[len(self.gaussian_collections)-1].pack_forget()
            del self.gaussian_collections[len(self.gaussian_collections)-1]

        while int(self.num_gaussian_collections.get()) > len(self.gaussian_collections):
            print "adding for final length of",int(self.num_gaussian_collections.get()) 
            self.gaussian_collections.append(GaussianField(self.gaussian_group.interior(), 
                                        tag_text = 'Gaussian mixture collection '+len(self.gaussian_collections).__repr__()))
            self.gaussian_collections[len(self.gaussian_collections)-1].create(name='gm'+len(self.gaussian_collections).__repr__(),
                                                                               callback=self.makeInitialParameters_preview)
            self.gaussian_collections[len(self.gaussian_collections)-1].pack(fill=Tkinter.X)

        self.makeInitialParameters_preview()


    def save(self):
        fname=self.parameters['initial_parameters_output_file']

        print "saving to",fname

        try:
            fd=open(fname,mode='w')
        except:
            print "Error opening file",fname
            return

        for l in self.text:
            fd.write(l+"\n")
        fd.close()
        print "saved"
        
    def makeInitialParameters_preview(self):
        #print "preview:"
        header_text=[]
        dcpt_text=[]
        dpmf_text=[]
        means_text=[]
        vars_text=[]
        gaussians_text=[]
        gm_text=[]

        # to do - make these user-definable
        mean_cell=' 0'
        var_cell=' 10'

        header_text.append('%gmtkPT: header begins')

        dcpt_text.append('% dense conditional probability tables')
        dcpt_text.append(len(self.dcpts).__repr__()+' % number of DCPTs')
        header_text.append('%gmtkPT: len(dcpts)='+len(self.dcpts).__repr__())
        n=0
        for dcpt in self.dcpts:
            dcpt_text.append(n.__repr__()+' % DCPT '+(n+1).__repr__()+' of '+len(self.dcpts).__repr__())
            dcpt_text.append(dcpt.name.get()+' % name')
            header_text.append('%gmtkPT: dcpt.name='+dcpt.name.get())
            header_text.append('%gmtkPT: dcpt.num_parents='+dcpt.num_parents.get())

            cards=''
            num_rows=1
            for p in dcpt.parents:
                cards=cards+p.get()+' '
                num_rows=num_rows*int(p.get())
                header_text.append('%gmtkPT: dcpt.parent.card='+p.get())

            header_text.append('%gmtkPT: dcpt.self_card='+dcpt.self_card.get())

            cards=cards+dcpt.self_card.get()
            dcpt_text.append(dcpt.num_parents.get()+' % num parents')
            dcpt_text.append(cards+' % parent cardinality(ies), self cardinality')
            row_text=""
            cell_text="%(val)4f " % {'val' : 1.0/float(dcpt.self_card.get())}
            for column in range(int(dcpt.self_card.get())):
                row_text=row_text+cell_text
            for row in range(num_rows):
                dcpt_text.append(row_text)

            n=n+1
            

        total_gc=0
        for gc in self.gaussian_collections:
            total_gc = total_gc + int(gc.num_gaussians.get())

        gm_text.append('% Gaussian mixture distributions')
        gm_text.append(total_gc.__repr__()+' % number of Gaussian mixture distributions')
        n=0
        for gc in self.gaussian_collections: # each gc is a Gaussian collection
            tn=0
            for g in range(int(gc.num_gaussians.get())): # each g is a Gaussian mixture distribution
                this_line=n.__repr__()+' '+gc.dimension.get() \
                           +' gmix_'+gc.name.get()+'_'+tn.__repr__() \
                           +' '+gc.num_mixes.get()+' ' \
                           +'dpmf_'+gc.name.get()+'_'+tn.__repr__()

                # make the dpmf needed for this mixture distribution
                dpmf_line=len(dpmf_text).__repr__()+' dpmf_'+gc.name.get()+'_'+tn.__repr__() \
                           +' '+gc.num_mixes.get()+' '
                cell_text="%(val)4f" % {'val' : 1.0/float(gc.num_mixes.get())}
                for gn in range(int(gc.num_mixes.get())):
                    dpmf_line=dpmf_line+' '+cell_text
                dpmf_text.append(dpmf_line)

                gnn=0
                for gn in range(int(gc.num_mixes.get())): # each gn is a single Gaussian
                    this_line=this_line+' gaussian_'+gc.name.get()+'_'+tn.__repr__()+'_'+gnn.__repr__()

                    # and now make this Gaussian, along with its mean and var
                    # the ' 0' is the type (0=diagonal covar, I think?)
                    gaussian_line=len(gaussians_text).__repr__() + ' ' + int(gc.dimension.get()).__repr__() \
                                   + ' 0' \
                                   +' gaussian_'+gc.name.get()+'_'+tn.__repr__()+'_'+gnn.__repr__() \
                                   +' mean_'+ gc.name.get()+'_'+tn.__repr__()+'_'+gnn.__repr__() \
                                   +' var_' + gc.name.get()+'_'+tn.__repr__()+'_'+gnn.__repr__()

                    mean_line=len(means_text).__repr__() + \
                               ' mean_'+ gc.name.get()+'_'+tn.__repr__()+'_'+gnn.__repr__()+' '+int(gc.dimension.get()).__repr__()
                    var_line=len(vars_text).__repr__() + \
                              ' var_'+ gc.name.get()+'_'+tn.__repr__()+'_'+gnn.__repr__()+' '+int(gc.dimension.get()).__repr__()
                    for d in range(int(gc.dimension.get())):
                        mean_line=mean_line+mean_cell
                        var_line=var_line+var_cell

                    gaussians_text.append(gaussian_line)
                    means_text.append(mean_line)
                    vars_text.append(var_line)
                    gnn=gnn+1
                                   

                gm_text.append(this_line)

                n=n+1
                tn=tn+1


        n=len(gaussians_text)
        prepend(gaussians_text,n.__repr__()+' % number of individual Gaussian components')
        prepend(gaussians_text,'% Gaussian components')

        n=len(dpmf_text)
        prepend(dpmf_text,n.__repr__()+'% number of DPMFs')
        prepend(dpmf_text,'% dense probability mass functions')

        n=len(means_text)
        prepend(means_text,n.__repr__())
        prepend(means_text,'% Means')

        n=len(vars_text)
        prepend(vars_text,n.__repr__())
        prepend(vars_text,'% Variances')

        self.text=[]
        self.text.append('% created by gmtkPT')
        self.text.append('')

        self.text.append('% ------------------------------')
        header_text.append('%gmtkPT: header ends')
        self.text.extend(header_text)

        self.text.append('% ------------------------------')
        self.text.extend(dpmf_text)
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.append('% sparse probability mass functions')
        self.text.append('0')
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.extend(means_text)
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.extend(vars_text)
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.append('% dlink matrices')
        self.text.append('0')
        self.text.append('% weight matrices')
        self.text.append('0')
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.extend(dcpt_text)
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.extend(gaussians_text)
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.extend(gm_text)
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.append('% switching mixture of Gaussians')
        self.text.append('0')
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.append('% logistic regression switching mixture of Gaussians')
        self.text.append('0')
        self.text.append('')

        self.text.append('% ------------------------------')
        self.text.append('% MLP-based switching mixture of Gaussians')
        self.text.append('0')

        # send text to preview field - TO DO: preserve scroll position
        #lo,hi=self.preview_box_scrollbar.get()
        #print "saving",lo,hi
        y=self.preview_box.yview()
        ln=1
        for l in self.text:
            self.preview_box.delete(ln.__repr__()+'.0',ln.__repr__()+'.end') # delete a line
            self.preview_box.insert(ln.__repr__()+'.0',l+"\n")
            ln=ln+1
        # delete any remaining lines
        self.preview_box.delete(ln.__repr__()+'.end',Tkinter.END)
            
        #self.preview_box_scrollbar.set(lo,hi)
        self.preview_box.yview(Tkinter.MOVETO,y[0])

            
    def clean_up(self):
        print "Cleaning up on exit..."
        self.kill_local()
        #self.kill_remote()

    def start_timer(self):
        # to fix:
        # check whether the agenda has self-terminated
        die_flag=True
        if self.agendas:
            for agenda in self.agendas:
                if self.agendas[agenda] and self.agendas[agenda].terminate_now != True:
                    die_flag=False
            if die_flag:
                print "All agendas terminated, dying"
                self.die()
            self.update()
            self.timer_id=self.update_button.after(self.interval*1000, self.start_timer)


    def update(self):
        if not self.gui_initialised:
            return

        cur_tab=self.notebook.getcurselection()

        if self.re_agenda.match(cur_tab):
            self.update_agenda(cur_tab)
        elif cur_tab == "Monitor":
            self.update_monitor()
        #else:
        #    print "Unknown tab in focus!"

    def update_agenda(self,name):
        print ">>>>>>>>>>>>> starting update_agenda"
        # from tab name, work out which agenda is being displayed on this tab
        agenda=self.agendas[name]

        #agenda.printme_brief()

        #print "update_agenda: on tab",name," we have agenda",agenda.getName()

        num_total=0
        num_completed=0
    
        # replace items, rather than emptying/refilling list, to maintain sort order and scroll position
        # this is too slow - needs speeding up and maybe threading ?
        if self.lists[name]:
            #cur_list0=list(self.lists[name].lists[0].get(0,Tkinter.END))
            
            if False: # cur_list0 != []:
                print ">>>>>>>>>>>>>>>>>> getting brief summary"

                (num_completed,num_total)=agenda.return_brief_summary(order=cur_list0,
                                                                      list=self.lists[name],
                                                                      num_total=num_total,
                                                                      num_completed=num_completed)
                print ">>>>>>>>>>>>>>>>>> got brief summary"


                if False:
                    # make a hash table of position in current list
                    pos_in_list={}
                    index=0
                    for i in cur_list0:
                        pos_in_list[i]=index
                        index=index+1

                    # update current list using new summary, replacing
                    # existing jobs at same position in list, inserting
                    # new jobs
                    for i in summary:
                        num_total=num_total+1
                        if i[1] == 'C':
                            num_completed=num_completed+1
                            try:
                                pos=pos_in_list[i[0]] #cur_list0.index(i[0]) # list may be not be sorted by key
                                #if self.lists[name].get(pos,pos) != [i]:
                                self.lists[name].delete(pos,pos)
                                self.lists[name].insert(pos,i)
                        
                            except ValueError:
                                self.lists[name].insert(Tkinter.END,i)

            else: # empty list - just insert all items
                (num_completed,num_total)=agenda.return_brief_summary(list=self.lists[name],
                                                                      num_total=num_total,
                                                                      num_completed=num_completed)
                #for i in agenda.return_brief_summary(self.lists[name]):
                #    num_total=num_total+1
                #    if i[1] == 'C':
                #        num_completed=num_completed+1
                #    self.lists[name].insert(Tkinter.END,i)

        if self.bestboxes[name]:
            self.bestboxes[name].delete('1.0',Tkinter.END)
            self.bestboxes[name].insert(Tkinter.END,agenda.return_current_best())
        print "<<<<<<<<<<<<< finishing update_agenda"

        # now update the progress bar
        self.progress_bars[name].updateProgress(num_completed,num_total)


        # and the estimated time to go
        if num_completed>0 and num_total>0 and agenda.start_time != None:

            this_time=time.time()
            so_far=this_time-agenda.start_time
            end_time=agenda.start_time + so_far*float(num_total)/float(num_completed)

            self.remaining_time[name].set('Estimated finish '
                                          +time.asctime(time.localtime(end_time)))

            if agenda.end_time != None:
                if end_time > agenda.end_time:
                    self.remaining_time[name].set('Estimated finish '
                                                  +time.asctime(time.localtime(end_time))
                                                  +' but will be aborted at '
                                                  +time.asctime(time.localtime(agenda.end_time)))

            

        else:
            self.remaining_time[name].set('Remaining time cannot be estimated')

        # simplest to flash all buttons
        for b in self.agenda_update_buttons:
            b.flash()

        print "Progress:",num_completed,num_total


    def update_monitor(self):
        self.update_local()
        if self.gui_initialised:
            self.monitor_update_button.flash()

        # cctrl can take many seconds to return, so thread this (to do??)
        self.remote_group_label_text.set('cctrl is running....')
        self.update_remote()
        self.remote_group_label_text.set('From cctrl')


        
    def update_local(self):

        fd = os.popen('ps -el')
        
        header_line=fd.readline()
        ps=fd.readlines()
        fd.close()

        # find which column the UID is in, etc
        s=header_line.split()
        uid_col=s.index("UID")
        pid_col=s.index("PID")
        ppid_col=s.index("PPID")

        lines=[]
        self.local_pids=[]
        if (not self.monitor_page) \
               or (self.local_process_monitor_users == None) \
               or (self.local_process_monitor_users.get() == "children"):

            # find also grandchildren, great-grandchildren,.....
            num_children=-1
            while num_children != len(self.local_pids):
                num_children = len(self.local_pids)
                for line in ps:
                    s=line.split()
                    if (int(s[ppid_col]) == int(self.this_pid)) or (self.local_pids.count(int(s[ppid_col])) > 0):
                        if self.local_pids.count(int(s[pid_col])) == 0:
                            self.local_pids.append(int(s[pid_col]))
                            lines.append(line)

        elif self.monitor_page and (self.local_process_monitor_users.get() == "self"):
            for line in ps:
                s=line.split()
                if int(s[uid_col]) == int(self.this_uid):
                    #self.local_pids.append(s[pid_col]) -- too dangerous: don't want the option of killing ALL these!!
                    lines.append(line)
        elif self.monitor_page and (self.local_process_monitor_users.get() == "all"):
            for line in ps:
                #self.local_pids.append(s[pid_col])
                lines.append(line)
        else:
            print "unknown local_process_monitor_users value:",self.local_process_monitor_users.get()
            
        if self.gui_initialised and self.monitor_page:
            self.localbox.delete('1.0',Tkinter.END)
            self.localbox.insert(Tkinter.END,header_line)
            for l in lines:
                self.localbox.insert(Tkinter.END,l)



    def kill_local(self):
            
        self.update_local()
        if len(self.local_pids) > 0:
            for p in self.local_pids:
                try:
                    print "killing",int(p)
                    os.kill(int(p),signal.SIGKILL)
                except OSError:
                    print "didn't kill",p
        else:
            print "Couldn't find any children"

        self.update_local()

    def kill_remote(self):

        # not good enough - hangs due to Catch-22 (cannot rexport to busy machines....)
            
        self.update_remote()
        if len(self.remote_jobs) > 0:
            for p in self.remote_jobs:
                try:
                    print "killing",p[1],"on node",p[0]
                    fd=os.popen('rexport -attr '+p[0]+' kill -9 '+p[1])
                    fd.close()
                except OSError:
                    print "didn't kill",p
        else:
            print "Couldn't find any remote processes"

        self.update_remote()

    def update_remote(self):

        fd = os.popen('/usr/local/sbin/cctrl -jobs -all')
        header_line=fd.readline()
        header_line=fd.readline()
        header_line=fd.readline() # third line has headers
        ps=fd.readlines()
        fd.close()

        # find which column the username is in, etc
        s=header_line.split()
        user_col=s.index("USER")
        host_col=s.index("HOST")
        pid_col=s.index("PID")

        lines=[]
        self.remote_jobs=[]
        if self.remote_process_monitor_users.get() == "self":
            for line in ps:
                s=line.split()
                if s[user_col] == self.this_username:
                    lines.append(line)
                    self.remote_jobs.append([s[host_col],s[pid_col]])
        elif self.remote_process_monitor_users.get() == "all":
            for line in ps:
                lines.append(line)
        else:
            print "unknown remote_process_monitor_users value!"
            
        if self.gui_initialised:
            self.remotebox.delete('1.0',Tkinter.END)
            self.remotebox.insert(Tkinter.END,header_line)
            for l in lines:
                self.remotebox.insert(Tkinter.END,l)

        print "here are the remote jobs",self.remote_jobs

    def select_one(self,e):
        # work out which tab we are
        name=self.notebook.getcurselection()
        if len(self.lists[name].curselection())>0:
            cs=[ int(x) for x in self.lists[name].curselection() ]

            # now grab first column to get job key
            jkey=int(self.lists[name].get(cs[0])[0])
            #print "jkey is",jkey
            t=self.agendas[name].jobs[jkey].output
            #print "got",t
            self.output_group_label_texts[name].set('Job '+jkey.__repr__()+' (double click job to update)')
            self.show_output(output=t,command=self.agendas[name].jobs[jkey].command)
            
    def show_output(self,command="",output=[]):
        # work out which tab we are
        name=self.notebook.getcurselection()
        self.cmdboxes[name].delete('1.0',Tkinter.END)
        # make command look a bit nicer
        self.cmdboxes[name].insert(Tkinter.END,command.replace(" -","\n -"))
        
        self.textboxes[name].delete('1.0',Tkinter.END)
        for l in output:
            self.textboxes[name].insert(Tkinter.END,l)
            
    def die(self):

        # grab the lock from the agenda to stop it exporting any more jobs
        print "GUI is dying, but waiting for agenda locks"
        # to fix
        if self.agendas:
            for agenda in self.agendas:
                if self.agendas[agenda] and self.agendas[agenda].thread_lock:
                    self.agendas[agenda].thread_lock.acquire()
        
        self.clean_up()

        if self.agendas:
            for agenda in self.agendas:
                if self.agendas[agenda]:
                    self.agendas[agenda].terminate_now=True
                    if self.agendas[agenda].thread_lock:
                        self.agendas[agenda].thread_lock.release()
        print "GUI ready to exit"
        
    def quit_button_pressed(self):
        print "Quit button pressed"
        self.die()
        if self.master:
            self.master.quit()
        print "GUI exiting now"

## ---------------------------------------------------------------

class DcptField(Pmw.Group):

    saved_callback=None

    def create(self, name='Name', callback=None):
        self.saved_callback=callback
        print "creating a DcptField with callback",callback
        self.num_parents=None
        self.name=Pmw.EntryField(self.interior(),
                                 command=callback,
                                 labelpos = 'n',
                                 label_text = 'name',
                                 validate = {'validator' : None,
                                             'min' : 0, 'max' : 99,},
                                 value=name)
        self.name.grid(column=0,row=0)


        self.num_parents = myPmwCounter(self.interior(),
                                 callback=self.make_parents,
                                 datatype='integer',
                                 labelpos = 'n',
                                 label_text = '# parents',
                                 orient = 'horizontal',
                                 entry_width = 4,
                                 entryfield_command=self.make_parents,
                                 entryfield_value = 1,
                                 entryfield_validate = {'validator' : 'integer',
                                                        'min' : 0, 'max' : 9999})
        self.num_parents.grid(column=1,row=0)

        self.self_card = myPmwCounter(self.interior(),
                                 callback=callback,
                                 datatype='integer',
                                 labelpos = 'n',
                                 label_text = 'self card',
                                 orient = 'horizontal',
                                 entry_width = 4,
                                 entryfield_command=callback,
                                 entryfield_value = 1,
                                 entryfield_validate = {'validator' : 'integer',
                                                        'min' : 1, 'max' : 9999})
        self.self_card.grid(column=2,row=0)

        
        self.parents=[]
        self.make_parents(callback)

    def make_parents(self,callback=None):
        if not callback:
            callback=self.saved_callback
        if not self.num_parents:
            print "no num_parents"
            return
        
        while int(self.num_parents.get()) < len(self.parents):
            print "forgetting down to",int(self.num_parents.get())
            self.parents[len(self.parents)-1].grid_forget()
            del self.parents[len(self.parents)-1]

        while int(self.num_parents.get()) > len(self.parents):
            print "adding for final length of",int(self.num_parents.get()) 
            self.parents.append(myPmwCounter(self.interior(),
                                    callback=callback,
                                    datatype='integer',
                                    labelpos = 'n',
                                    label_text = 'p'+len(self.parents).__repr__()+' card',
                                    orient = 'horizontal',
                                    entry_width = 4,
                                    entryfield_value = 1,
                                    entryfield_command=callback,
                                    entryfield_validate = {'validator' : 'integer',
                                                           'min' : 1, 'max' : 9999}))

            self.parents[len(self.parents)-1].grid(column=len(self.parents)+2,row=0)

        apply(callback)





class GaussianField(Pmw.Group):

    def create(self,name='Name', callback=None):
        self.name=Pmw.EntryField(self.interior(), 
                                 command=callback,
                                 validate = {'validator' : None,
                                             'min' : 0, 'max' : 99,},
                                 value=name)
        self.name.grid(column=0,row=0)

        self.num_gaussians = myPmwCounter(self.interior(),
                                          callback=callback,
                                          datatype='integer',
                                          labelpos = 'n',
                                          label_text = ' # GMMs',
                                          orient = 'horizontal',
                                          entry_width = 4,
                                          entryfield_value = 1,
                                          entryfield_command=callback,
                                          entryfield_validate = {'validator' : 'integer',
                                                            'min' : 1, 'max' : 9999})
        self.num_gaussians.grid(column=1,row=0)

        self.dimension = myPmwCounter(self.interior(),
                                      callback=callback,
                                      datatype='integer',
                                      labelpos = 'n',
                                      label_text = ' obs dim',
                                      orient = 'horizontal',
                                      entry_width = 3,
                                      entryfield_value = 1,
                                      entryfield_command=callback,
                                      entryfield_validate = {'validator' : 'integer',
                                                            'min' : 1, 'max' : 999})
        self.dimension.grid(column=2,row=0)
        
        self.num_mixes = myPmwCounter(self.interior(),
                                       callback=callback,
                                       datatype='integer',
                                       labelpos = 'n',
                                       label_text = ' # mix comps / GMM',
                                       orient = 'horizontal',
                                       entry_width = 2,
                                       entryfield_value = 1,
                                       entryfield_command=callback,
                                       entryfield_validate = {'validator' : 'integer',
                                                              'min' : 1, 'max' : 99})
        self.num_mixes.grid(column=3,row=0)





class ProgressBar:

    "This comes from John Grayson's book Python and Tkinter programming"

    def __init__(self, master=None, orientation="horizontal",
                 min=0, max=100, width=100, height=18,
                 doLabel=1, appearance="sunken",
                 fillColor="blue", background="gray",
                 labelColor="yellow", labelFont="Verdana",
                 labelText="", labelFormat="%d%%",
                 value=50, bd=2):
        # preserve various values
        self.master=master
        self.orientation=orientation
        self.min=min
        self.max=max
        self.width=width
        self.height=height
        self.doLabel=doLabel
        self.fillColor=fillColor
        self.labelFont= labelFont
        self.labelColor=labelColor
        self.background=background
        self.labelText=labelText
        self.labelFormat=labelFormat
        self.value=value
        self.frame=Tkinter.Frame(master, relief=appearance, bd=bd)
        self.canvas=Tkinter.Canvas(self.frame, height=height, width=width, bd=0,
                           highlightthickness=0, background=background)
        self.scale=self.canvas.create_rectangle(0, 0, width, height,
                                                fill=fillColor)
        self.label=self.canvas.create_text(self.canvas.winfo_reqwidth()
/ 2,
                                           height / 2, text=labelText,
                                           anchor="c", fill=labelColor,
                                           font=self.labelFont)
        self.update()
        self.canvas.pack(side='top', fill='x', expand='no')

    def updateProgress(self, newValue, newMax=None):
        if newMax:
            self.max = newMax
        self.value = newValue
        self.update()

    def update(self):
        # Trim the values to be between min and max
        value=self.value
        if value > self.max:
            value = self.max
        if value < self.min:
            value = self.min
        # Adjust the rectangle
        if self.orientation == "horizontal":
            self.canvas.coords(self.scale, 0, 0,
              float(value) / self.max * self.width, self.height)
        else:
            self.canvas.coords(self.scale, 0,
                               self.height - (float(value) / 
                                              self.max*self.height),
                               self.width, self.height)
        # Now update the colors
        self.canvas.itemconfig(self.scale, fill=self.fillColor)
        self.canvas.itemconfig(self.label, fill=self.labelColor)
        # And update the label
        if self.doLabel:
            if value:
                if value >= 0:
                    pvalue = int((float(value) / float(self.max)) * 
                                   100.0)
                else:
                    pvalue = 0
                self.canvas.itemconfig(self.label, text=self.labelFormat
                                         % pvalue)
            else:
                self.canvas.itemconfig(self.label, text='')
        else:
            self.canvas.itemconfig(self.label, text=self.labelFormat %
                                   self.labelText)
        self.canvas.update_idletasks()


