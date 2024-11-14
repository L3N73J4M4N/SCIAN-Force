from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
from forces import Forces3D
import numpy as np
import pandas as pd
from PIL import ImageTk
from icons import icon_force_karina


class Inter3D:
    def __init__(self):
        self.obj = None  # objeto con el que se trabajará
        win = Tk()  # ventana
        x_win, y_win = 370, 465
        width = x_win - 20
        win.minsize(x_win, y_win)  # tamaño de la ventana
        win.title("SCIAN-Force (v2024.11.08)")  # título de la ventana
        icon16, icon32 = icon_force_karina()
        icon16, icon32 = ImageTk.PhotoImage(icon16), ImageTk.PhotoImage(icon32)
        win.iconphoto(False, icon32, icon16)
        x_screen, y_screen = win.winfo_screenwidth(), win.winfo_screenheight()
        x = x_screen // 2 - x_win // 2
        y = y_screen // 2 - y_win // 2
        win.geometry(str(x_win) + 'x' + str(y_win) + '+' + str(x) + '+' + str(y))
        win.resizable(None, None)
        self.use_stress = BooleanVar(value=True)  # mostrar estrés o curvaturas

        def entry_data():
            value1.delete(0, 'end')
            value2.delete(0, 'end')
            value3.delete(0, 'end')
            value4.delete(0, 'end')
            mean = np.mean(self.obj.meanCurvature)
            std = np.std(self.obj.meanCurvature)
            maxi = np.max(self.obj.meanCurvature)
            mini = np.min(self.obj.meanCurvature)
            value1.insert(0, str(np.round(mean, 10)))
            value2.insert(0, str(np.round(std, 10)))
            value3.insert(0, str(np.round(maxi, 10)))
            value4.insert(0, str(np.round(mini, 10)))

        def openfile():
            filepath = filedialog.askopenfilename(initialdir='Desktop',
                                                  title='Select .OFF / .OBJ file',
                                                  filetypes=(('object files', '*.obj *.off'),
                                                             ('obj files', '*.obj'),
                                                             ('off files', '*.off'),
                                                             ('all files', '*.*')))
            if len(filepath) > 0:
                try:
                    self.obj = Forces3D(filepath, interval=interval.get(), rotate=True)
                except ValueError:
                    self.obj = Forces3D(filepath, rotate=True)
                self.obj.fix_mesh()
                string = interval.get().split('x')
                string = '(' + string[0] + ' x ' + string[1] + ' x ' + string[2] + ')'
                label_down.configure(text=string)
                name = self.obj.name
                if len(self.obj.name) > 30:
                    name = self.obj.name[0:30] + '\n' + self.obj.name[29:]
                if len(self.obj.name) > 60:
                    name = self.obj.name[0:30] + '\n' + self.obj.name[29:60] + '\n' + self.obj.name[59:]
                obj_name.configure(text=name)
                radius = str(self.obj.get_radius())
                volume = str(self.obj.get_volume())
                if len(radius) >= 10:
                    radius = radius[0:10]
                if len(volume) >= 10:
                    volume = volume[0:10]
                value5.delete(0, 'end')
                value6.delete(0, 'end')
                value5.insert(0, radius)
                value6.insert(0, volume)

        def laplace():
            self.obj.optimized_laplacian()
            try:
                max_value = float(maxvalue.get())
            except ValueError:
                max_value = None
            try:
                min_value = float(minvalue.get())
            except ValueError:
                min_value = None
            self.gamma = 1 if gamma.get() == '' else float(gamma.get())
            if self.use_stress.get():
                self.obj.stress(self.gamma)
                print('OK')
            self.obj.limit_curvatures(min_value, max_value)
            self.obj.plot_go(max_value=max_value)
            entry_data()

        def rus():
            self.obj.rusinkiewicz_curvature()
            try:
                max_value = float(maxvalue.get())
            except ValueError:
                max_value = None
            try:
                min_value = float(minvalue.get())
            except ValueError:
                min_value = None
            self.gamma = 1 if gamma.get() == '' else float(gamma.get())
            if self.use_stress.get():
                self.obj.stress(self.gamma)
            self.obj.limit_curvatures(min_value, max_value)
            self.obj.plot_go(max_value=max_value)
            entry_data()

        # object
        select_obj = LabelFrame(win, text='Object')
        select_obj.pack(pady=2)
        button_select = Button(select_obj, text='Select object file', command=openfile)
        obj_name = Label(select_obj, text='')
        interval = Entry(select_obj)
        interval.insert(0, '1x1x1')
        label_interval = Label(select_obj, text='Voxel size [μm³]:')
        label_down = Label(select_obj, text='(width x height x depth)')
        label_gamma = Label(select_obj, text='γ [mN/m]:')
        gamma = Entry(select_obj)
        use_stress = Checkbutton(select_obj, text='Show stress', variable=self.use_stress)

        button_select.grid(row=2, columnspan=2)
        obj_name.grid(row=3, columnspan=2)
        interval.grid(row=0, column=1, padx=10)
        label_interval.grid(row=0, column=0)
        label_down.grid(row=1, columnspan=2)
        label_gamma.grid(row=4, column=0)
        gamma.grid(row=4, column=1)
        use_stress.grid(row=5, columnspan=2)

        # colorbar
        cb = LabelFrame(win, text='Clipping', width=width)
        cb.pack(pady=2)
        maxvalue_label = Label(cb, text='Max value:')
        minvalue_label = Label(cb, text='Min value:')
        maxvalue = Entry(cb, width=15)
        minvalue = Entry(cb, width=15)
        space = Label(cb, text='')
        space.grid(column=3, row=0, padx=9)
        space = Label(cb, text='')
        space.grid(column=0, row=0, padx=9)
        maxvalue_label.grid(column=1, row=0)
        minvalue_label.grid(column=1, row=1)
        maxvalue.grid(column=2, row=0, padx=20)
        minvalue.grid(column=2, row=1)

        # plot
        plot = LabelFrame(win, text='Plot', width=width)
        plot.pack(pady=2)
        button_laplace = Button(plot, text='Laplacian method', command=laplace)
        button_rus = Button(plot, text='Rusinkiewicz method', command=rus)

        # pack de los botones
        space1 = Label(plot)
        space2 = Label(plot)
        button_laplace.grid(column=1, row=0)
        button_rus.grid(column=1, row=1)
        space1.grid(column=0, row=0, padx=28)
        space2.grid(column=2, row=0, padx=28)

        # valores
        res = LabelFrame(win, text='Results', width=width)
        res5 = Label(res, text='radius [μm]:')
        value5 = Entry(res)
        res6 = Label(res, text='volume [μm³]:')
        value6 = Entry(res)
        res1 = Label(res, text=r'mean H/σ :')
        value1 = Entry(res)
        res2 = Label(res, text='std H/σ :')
        value2 = Entry(res)
        res3 = Label(res, text='max H/σ :')
        value3 = Entry(res)
        res4 = Label(res, text='min H/σ :')
        value4 = Entry(res)

        # guardar los resultados de las curvaturas
        def save():
            str_save = 'normal_stress' if self.use_stress else 'mean_curvature'
            dicty = {'vertex': [], 'coord': [], str_save: []}
            for i in range(0, len(self.obj.meanCurvature)):
                dicty['vertex'].append(i)
                coord = self.obj.mesh[0][i]
                dicty['coord'].append((coord[0], coord[1], coord[2]))
                dicty[str_save].append(self.obj.meanCurvature[i])
            df = pd.DataFrame(dicty)
            name_path = filedialog.asksaveasfilename(initialdir="Desktop",
                                                     initialfile=self.obj.name[0:-4] + '-' + self.obj.method[19:],
                                                     defaultextension='.xlsx',
                                                     title="Save as",
                                                     filetypes=(("xlsx files", "*.xlsx"),
                                                                ("all files", "*.*")))
            if len(name_path) > 0:
                if name_path[-5::] != '.xlsx':
                    name_path += '.xlsx'
                df.to_excel(name_path, sheet_name='results', index=False)

        save_curvatures = Button(res, text='Save stresses or curvatures', command=save)

        res.pack(pady=2)
        res5.grid(column=1, row=4)
        res6.grid(column=1, row=5)
        value5.grid(column=2, row=4)
        value6.grid(column=2, row=5)
        res1.grid(column=1, row=0)
        res2.grid(column=1, row=1)
        res3.grid(column=1, row=2)
        res4.grid(column=1, row=3)

        value1.grid(column=2, row=0)
        value2.grid(column=2, row=1)
        value3.grid(column=2, row=2)
        value4.grid(column=2, row=3)
        save_curvatures.grid(row=6, columnspan=4)

        space3 = Label(res)
        space4 = Label(res)
        space3.grid(column=0, row=0, padx=7)
        space4.grid(column=3, row=0, padx=7)

        # iniciar ventana
        win.mainloop()


Inter3D()
