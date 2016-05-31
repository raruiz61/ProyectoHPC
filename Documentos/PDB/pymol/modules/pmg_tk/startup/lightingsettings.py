'''
Lighting Settings Plugin

(c) 2013 Schrodinger Inc.
'''

import Tkinter
from pymol import cmd, plugins

def __init_plugin__(self=None):
    plugins.addmenuitem('Lighting Settings', lightingsettings)

def sliderentry(master, var, **kwargs):
    frame = Tkinter.Frame(master)
    scale = Tkinter.Scale(frame, showvalue=0, orient=Tkinter.HORIZONTAL,
            variable=var, **kwargs)
    entry = Tkinter.Entry(frame, textvariable=var, width=5)
    entry.pack(side=Tkinter.RIGHT)
    scale.pack(side=Tkinter.LEFT, expand=1, fill=Tkinter.X)
    return frame

def preset_default():
    cmd.set('ambient', 0.14)
    cmd.set('direct', 0.45)
    cmd.set('spec_direct', 0)
    cmd.set('spec_direct_power', 55.)
    cmd.set('light_count', 2)
    cmd.set('shininess', 55.)
    cmd.set('reflect', 0.45)
    cmd.set('spec_count', -1)
    cmd.set('spec_power', -1.)
    cmd.set('spec_reflect', -1.)
    cmd.set('specular', 1)
    cmd.set('specular_intensity', 0.5)
    cmd.set('ambient_occlusion_mode', 0)
    cmd.set('ambient_occlusion_scale', 25.0)
    cmd.set('ambient_occlusion_smooth', 10)
    cmd.set('power', 1.)
    cmd.set('reflect_power', 1.)

def preset_metal():
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.) # diffuse
    cmd.set('spec_direct', 0)
    cmd.set('shininess', 51.) # same as spec_power
    cmd.set('reflect', 0.5) # diffuse
    cmd.set('spec_count', -1)
    cmd.set('spec_reflect', -1.)
    cmd.set('specular', 1)
    cmd.set('specular_intensity', 0.5)  # same as specular

def preset_plastic():
    cmd.set('ambient', 0.)
    cmd.set('direct', 0.2) # diffuse
    cmd.set('spec_direct', 0)
    cmd.set('shininess', 32.) # same as spec_power
    cmd.set('reflect', 0.55) # diffuse
    cmd.set('spec_count', -1)
    cmd.set('spec_reflect', -1.)
    cmd.set('specular', 1)
    cmd.set('specular_intensity', 0.5)  # same as specular

def preset_rubber():
    cmd.set('ambient', 0.05)
    cmd.set('direct', 0.2) # diffuse
    cmd.set('spec_direct', 0)
    cmd.set('shininess', 10.) # same as spec_power
    cmd.set('reflect', 0.5) # diffuse
    cmd.set('spec_count', -1)
    cmd.set('spec_reflect', -1.)
    cmd.set('specular', 1)
    cmd.set('specular_intensity', 0.5)  # same as specular

def lightingsettings():
    self = Tkinter.Toplevel(plugins.get_tk_root())
    self.title('Lighting Settings')

    setting = plugins.get_pmgapp().skin.setting

    sliders = [
        "Diffuse Reflection",
        ('ambient', 0, 1, None),
        ('reflect', 0, 1, None),

        "Direct Light from Front",
        ('direct (+reflect)', 0, 1, None), # diffuse, coupled with "reflect"
        ('spec_direct', 0, 1, None),
        ('spec_direct_power', 0, 100, 1),

        "Free placeable directed Lights",
        ('light_count', 1, 8, 1),
        ('edit_light', 1, 7, 1),

        "Specular Reflection",
        ('spec_count', -1, 8, 1),
        # ('spec_power', -1, 200, 1), # deprecated since v1.5
        ('shininess', 0, 100, None), # same as spec_power
        ('spec_reflect', -0.01, 1, None),
        ('specular', 0, 1, None),
        ('specular_intensity (=specular)', 0, 1, None), # same as specular

        "Ambient Occlusion (Surface only)",
        ('ambient_occlusion_mode', 0, 2, 1),
        ('ambient_occlusion_scale', 1.0, 50., None),
        ('ambient_occlusion_smooth', 1, 20, 1),

        "Ray trace only",
        ('power', -10, 100, 1),
        ('reflect_power', 0, 10, 1),
    ]

    frame = Tkinter.Frame(self)

    Tkinter.Label(frame, text="Presets:", fg="red").pack(side=Tkinter.LEFT)
    Tkinter.Button(frame, text="Default", command=preset_default).pack(side=Tkinter.LEFT)
    Tkinter.Button(frame, text="Metal", command=preset_metal).pack(side=Tkinter.LEFT)
    Tkinter.Button(frame, text="Plastik", command=preset_plastic).pack(side=Tkinter.LEFT)
    Tkinter.Button(frame, text="Rubber", command=preset_rubber).pack(side=Tkinter.LEFT)

    frame.grid(row=0, columnspan=2, sticky=Tkinter.W)

    for i, item in enumerate(sliders, 1):
        if isinstance(item, str):
            Tkinter.Label(self, text=item, fg="blue").grid(
                    row=i, columnspan=2, sticky=Tkinter.W)
            continue

        name, min, max, res = item
        if res is None:
            res = 0.01 if (max - min < 100) else 0.1

        var = getattr(setting, name.split()[0])

        slider = sliderentry(self, var, from_=min, to=max,
                resolution=res, length=200)

        tlabel = Tkinter.Label(self, text=name)

        tlabel.grid(row=i, column=0, sticky=Tkinter.W)
        slider.grid(row=i, column=1, sticky=Tkinter.E)

# vi:expandtab:smarttab
