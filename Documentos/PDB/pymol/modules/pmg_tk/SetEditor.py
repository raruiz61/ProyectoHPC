#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2006 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-* Kenneth Lind
#Z* -------------------------------------------------------------------

# Master editor for all pymol settings
# (includes filter)

from Tkinter import *
import tkMessageBox

class MiscEntry(Entry):
    def set(self, value):
        self.delete(0, END)
        self.insert(END, value)

def bind_mousewheel(scrolledframe):
    # Enable scrolling with mouse wheel on Pmw.ScrolledFrame
    delta = 0.1

    def scrollUp(event=None):
        scrolledframe.yview('scroll', -delta, 'pages')

    def scrollDown(event=None):
        scrolledframe.yview('scroll', delta, 'pages')

    def bind_rec(other):
        other.bind('<Button-4>', scrollUp)
        other.bind('<Button-5>', scrollDown)
        for child in other.winfo_children():
            bind_rec(child)
    bind_rec(scrolledframe)

class SetEditor:

    def __init__(self, app, parent=None):
        import Pmw

        self.pymol = app.pymol
        self.cmd = app.pymol.cmd
        
        #====================
        # get a list of settings and create dictionary with setting:value items
        name_list = self.pymol.setting.get_name_list()
        name_list.sort()
        self.rows = []

        #====================
        # create window
        top = parent
        if not top:
            top = Toplevel()
            top.title( "PyMOL Settings" )
            top.geometry( "+100+200" )
            top.protocol('WM_DELETE_WINDOW', self.onClose)

        pw = PanedWindow(top, orient=VERTICAL)
        pw.pack(fill=BOTH, expand=1)

        #====================
        # help frame
        sf = Pmw.ScrolledFrame(pw, hull_height=100, usehullsize=1)
        sf.pack(side=BOTTOM, fill=X, expand=0)
        f3l = Label(sf.interior(), text="", anchor=W, wraplength=500, justify='left')
        f3l.pack(fill=X, expand=1)
        self.help_label = f3l
        bind_mousewheel(sf)

        f4 = Frame(pw)
        pw.add(f4)
        pw.add(sf)

        #====================
        # create filter frame
        f2 = Frame( f4, bd=1, relief=SUNKEN )
        f2.pack( side=BOTTOM, fill=X, expand=0 )
        l = Label( f2, text="Filter:", anchor=E )
        l.pack(side=LEFT)
        self.filter = Entry( f2 )
        self.filter.bind( "<KeyRelease>", self.onFilter )
        self.filter.pack(side=LEFT, fill=X, expand=1)
        self.filter.focus_set()
        self.filter_help = IntVar()
        c = Checkbutton(f2, text="Help Text", variable=self.filter_help,
                command=self.onFilter)
        c.pack(side=RIGHT)

        #====================
        # scrolled table
        sf = Pmw.ScrolledFrame(f4, vscrollmode='static', horizflex='expand')
        sf.pack(side=TOP, fill=BOTH, expand=1)
        f1 = sf.interior()
        for i, lab in enumerate(name_list):
            l = Label(f1, text=lab, anchor=E)
            l.grid(row=i, column=0, sticky=E)

            v_type, v_list = self.cmd.get_setting_tuple(lab)
            if v_type == 1:
                e = Checkbutton(f1, variable=getattr(app.setting, lab))
                e.grid(row=i, column=1, sticky=W)
            else:
                e = MiscEntry(f1)
                e.grid(row=i, column=1, sticky=W+E, padx=7)
                callback = lambda event, i=i: self.onSet(i)
                e.bind("<Return>", callback)
                e.bind("<FocusOut>", callback)
                e.set(self.cmd.get(lab))
            callback = lambda event, i=i: self.onFocus(offset=i)
            e.bind("<Enter>", callback)
            e.bind("<FocusIn>", callback)

            row = (lab, l, e)
            self.rows.append(row)
        f1.columnconfigure(1, weight=1)
        bind_mousewheel(sf)

        self.offset = 0
        self.top = top

    def onFocus(self, event=None, offset=-1):
        from pymol.setting_help import setting_help_dict

        text = '-- no help available --'
        self.offset = offset

        if offset >= 0:
            try:
                lab = self.rows[offset][0]
                text, stype, default, level = setting_help_dict[lab]
                level = ['global', 'object', 'object-state', 'atom', 'atom-state'][int(level)]
                if stype:
                    if default:
                        text = '%s = %s {default: %s} (%s):\n\n%s' % (lab, stype, default, level, text)
                    else:
                        text = '%s = %s (%s):\n\n%s' % (lab, stype, level, text)
                else:
                    text = '%s (%s):\n\n%s' % (lab, level, text)
            except LookupError:
                pass

        self.help_label.config(text=text)

    def onSet(self, offset):
        """set the pymol setting to value in the entry widget where
        the enter key was pressed
        """
        lab, _, entry = self.rows[offset]
        try:
            val = entry.get()
        except AttributeError:
            return
        if val == self.cmd.get(lab):
            return
        try:
            self.cmd.set( lab, val, quiet=0, log=1 )
        except:
            tkMessageBox.showerror('Error', 'Failed to set value for "%s"' % (lab))
        entry.set(self.cmd.get(lab))

    def onClose(self):
        self.onSet(self.offset)
        self.top.destroy()

    def onFilter(self, *event):
        """get list of pymol settings that match filter"""
        from pymol.setting_help import setting_help_dict

        val = self.filter.get().strip().lower()
        filter_help = self.filter_help.get()
        for lab, l, e in self.rows:
            remove = val not in lab
            if remove and filter_help:
                try:
                    text = setting_help_dict[lab][0].lower()
                    remove = val not in text
                except LookupError:
                    pass
            if remove:
                l.grid_remove()
                e.grid_remove()
            else:
                l.grid()
                e.grid()

