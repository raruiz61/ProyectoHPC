import os.path
from pymol.wizard import Wizard
from pymol import cmd, controlling

class Find(Wizard):
    '''
    Find an item in the object menu panel and scroll it to the top
    '''

    prev_pattern = ''

    def __init__(self, _self=cmd):
        Wizard.__init__(self, _self)

        self.pattern = Find.prev_pattern
        self.hitnum = 0
        self.hitcount = 0

        self.find_now()

    def get_event_mask(self):
        return Wizard.event_mask_key

    def do_key(self,k,x,y,m):
        if (m,k) in [
                (2,6),    # CTRL+F
                (0,13),   # Return
            ]:
            self.hitnum += 1
        elif (m,k) in [
                (1,13),   # SHIFT+Return
            ]:
            self.hitnum -= 1
        elif k == 9: # Tab
            self.pattern = os.path.commonprefix([n
                    for n in self.cmd.get_names('public')
                    if n.startswith(self.pattern)])
        elif (m,k) in [
                (1,8),      # SHIFT+Backspace
                (1,127),
            ]:
            self.pattern = ""
            self.hitnum = 0
        elif k in [8, 127]: # Backspace
            self.pattern = self.pattern[:-1]
            self.hitnum = 0
        elif k == 27:   # ESC
            self.cancel()
            return 1
        elif 32 < k < 127:
            self.pattern += chr(k)
            self.hitnum = 0
        else:
            return 1
        self.find_now()
        return 1

    def find_now(self):
        self.hitcount = controlling.find(self.pattern, self.hitnum)
        if self.hitcount:
            self.hitnum %= self.hitcount
        self.cmd.refresh_wizard()

    def get_prompt(self):
        self.prompt = [
            r"\090Type to find object or selection in panel",
            r"<Enter>        next hit",
            r"<Shift+Enter>  previous hit",
            r"<Tab>          auto complete",
            r"<Esc>          cancel",
            r"",
            r"\559Find:\999 " + self.pattern + "_"
            r"",
        ]
        if self.pattern:
            self.prompt += [
                "",
                r"\900No hit"
                if self.hitcount == 0 else
                r"\990Found %d hit(s)" % self.hitcount,
            ]
        return self.prompt

    def get_panel(self):
        return [
            [ 1, 'Find in Panel', '' ],
            [ 2, 'Cancel', 'cmd.get_wizard().cancel()' ]
            ]

    def cancel(self):
        Find.prev_pattern = self.pattern
        self.cmd.set_wizard()
        self.cmd.refresh_wizard()

