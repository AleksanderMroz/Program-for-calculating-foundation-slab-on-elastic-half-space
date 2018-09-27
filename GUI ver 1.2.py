
from tkinter import messagebox
import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg


import tkinter as tk
from tkinter import ttk

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


import BackEND as BE

LARGE_FONT = ("Verdana", 12)

d = plt.figure()
a = d.add_subplot(111)


class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        lbl2 = tk.Label(self, text="Parametry \n płyty:", font=("Arial Bold", 16))
        lbl2.grid(column=0, row=1)

        ##########################
        lbl3 = tk.Label(self, text="Wymiar a[m] ",
                        font=("Arial", 12))
        lbl3.grid(column=0, row=2)

        txt3 = tk.Entry(self, width=10)
        txt3.grid(column=1, row=2)
        ##########################
        lbl4 = tk.Label(self, text="Wymiar b[m] ",
                        font=("Arial", 12))
        lbl4.grid(column=0, row=3)

        txt4 = tk.Entry(self, width=10)
        txt4.grid(column=1, row=3)
        ##########################
        lbl5 = tk.Label(self, text="Grubość h[m] ",
                        font=("Arial", 12))
        lbl5.grid(column=0, row=4)

        txt5 = tk.Entry(self, width=10)
        txt5.grid(column=1, row=4)
        ##########################
        lbl6 = tk.Label(self, text="Poisson v[-] ",
                        font=("Arial", 12))
        lbl6.grid(column=0, row=5)

        txt6 = tk.Entry(self, width=10)
        txt6.grid(column=1, row=5)
        ##########################
        lbl7 = tk.Label(self, text="E [MPa] ",
                        font=("Arial", 12))
        lbl7.grid(column=0, row=6)

        txt7 = tk.Entry(self, width=10)
        txt7.grid(column=1, row=6)
        ##########################

        ###########################
        # *************************#
        ###########################
        lab0 = tk.Label(self, text="Parametry \ngruntu:",
                        font=("Arial Bold", 16))
        lab0.grid(column=4, row=1)

        lab1 = tk.Label(self, text="E [MPa] ",
                        font=("Arial", 12))
        lab1.grid(column=4, row=2)
        tx1 = tk.Entry(self, width=10)
        tx1.grid(column=5, row=2)
        ###########################
        lab2 = tk.Label(self, text="Poisson v [-] ",
                        font=("Arial", 12))
        lab2.grid(column=4, row=3)
        tx2 = tk.Entry(self, width=10)
        tx2.grid(column=5, row=3)

        ###########################
        # *************************#
        ###########################

        l0 = tk.Label(self, text="Siatkowanie \nelementu:",
                      font=("Arial Bold", 16))
        l0.grid(column=4, row=5)
        #######################################
        l1 = tk.Label(self, text="Podział boku a:\nna m odcinkow",
                      font=("Arial", 12))
        l1.grid(column=4, row=6)

        spin1 = tk.Spinbox(self, from_=4, to=80, width=3)
        spin1.grid(column=5, row=6)
        #######################################
        l1 = tk.Label(self, text="Podział boku b:\nna n odcinkow",
                      font=("Arial", 12))
        l1.grid(column=4, row=7)

        spin2 = tk.Spinbox(self, from_=4, to=80, width=3)
        spin2.grid(column=5, row=7)

        l2 = tk.Label(self, text="(Zaleca się podział\n m*n<=400)",
                      font=("Arial", 12))
        l2.grid(column=4, row=8)
        #######################################

        Sc = tk.Label(self, text="Podaj sciezke \n zapisu pliku z wynikami\n"
                                 "np: d:\FolderX\Wyniki1.txt",
                      font=("Arial", 12))
        Sc.grid(column=0, row=10)

        ttt1 = tk.Entry(self, width=10)
        ttt1.grid(column=1, row=10)

        btn0 = tk.Button(self, text="OBLICZENIA",
                         command=lambda: controller.Oblicz(txt3.get(), txt4.get(), txt5.get(), txt6.get(),
                                                           txt7.get(),
                                                           tx1.get(), tx2.get(), spin1.get(), spin2.get(), controller.F,
                                                           ttt1.get(), controller.D))
        btn0.grid(column=6, row=2)

        l3 = tk.Label(self, text="Wprowadź obciążenie:",
                      font=("Arial", 12))
        l3.grid(column=6, row=3)
        btn = tk.Button(self, text="Obciążenia", command=lambda: controller.show_frame(PageObciążenia))
        btn.grid(column=6, row=4)

        l4 = tk.Label(self, text="Wyniki:",
                      font=("Arial", 12))
        l4.grid(column=6, row=6)
        btn1 = tk.Button(self, text="Siły kontaktowe", command=lambda: controller.show_frame(PageR))
        btn1.grid(column=6, row=7)

        btn2 = tk.Button(self, text="Osiadania", command=lambda: controller.show_frame(PageOsiadania))
        btn2.grid(column=7, row=7)

        btn3 = tk.Button(self, text="MXX", command=lambda: controller.show_frame(PageMXX))
        btn3.grid(column=6, row=8)

        btn4 = tk.Button(self, text="MXY", command=lambda: controller.show_frame(PageMXY))
        btn4.grid(column=7, row=8)

        btn5 = tk.Button(self, text="MYY", command=lambda: controller.show_frame(PageMYY))
        btn5.grid(column=6, row=9)


class SeaofBTCapp(tk.Tk):
    a0, b0, m, n = 0, 0, 0, 0
    Krotka = None
    F = [[0, 0, 0]]
    D = [[10000000000000000, 0, 0, 0]]
    MXX = "BLAH"
    MYY = "MEH"
    MXY = "WHAAAAT"
    R = "WUT"
    W = "NOPE"

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.wm_title(self, "Praca magisterska: Program do obliczania płyt prostokątnych na półprzestrzeni sprężystej")

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu2 = tk.Menu(menubar, tearoff=0)
        filemenu3 = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Zapisz", command=lambda: clicked3())
        filemenu.add_command(label="Wczytaj", command=lambda: clicked3())
        filemenu.add_separator()
        filemenu.add_command(label="Zakończ program", command=quit)

        filemenu2.add_command(label="O programie:", command=lambda: clicked1())
        filemenu2.add_command(label="O mnie:", command=lambda: clicked2())
        filemenu2.add_command(label="Poradnik:", command=lambda: clicked4())

        filemenu3.add_command(label="Wymiary", command=lambda: self.show_frame(HELPSTART))
        filemenu3.add_command(label="Obciążenia", command=lambda: self.show_frame(HELPOBC))

        menubar.add_cascade(label="Plik", menu=filemenu)
        menubar.add_cascade(label="O mnie", menu=filemenu2)
        menubar.add_cascade(label="Help", menu=filemenu3)

        tk.Tk.config(self, menu=menubar)

        self.frames = {}
        for F in (StartPage, PageObciążenia, PageMXX, PageMYY, PageOsiadania, PageMXY, PageR, HELPOBC, HELPSTART):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

    def Reboot(self, cont):
        frame = self.frames[cont]
        frame.destroy()
        self.frames = {}
        for F in (StartPage, PageObciążenia, PageMXX):
            frame = F(self.container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame(StartPage)

    def Calculate(self, *args):
        print(*args)
        print(self.MXX)
        self.MXX = KAL.ObliczCos(41, 41)
        print(self.MXX)

    def Oblicz(self, *args):
        self.a0 = args[0]
        self.b0 = args[1]
        self.m = args[7]
        self.n = args[8]
        self.Krotka = (self.a0, self.b0, self.m, self.n)
        print(self.Krotka, type(self.Krotka), self.Krotka[0], type(self.Krotka[0]))
        self.R, self.W, self.MXX, self.MXY, self.MYY = BE.Oblicz(*args)
        print("TESCIOR")
        print(np.shape(self.W))


def clicked1():
    messagebox.showinfo("O programie",
                        "Program komputerowy powstał w ramach pracy magisterskiej. Służy do wyznaczania naprezen kontaktowych  oraz osiadan plyty posadowionej na polprzestrzeni sprezystej. \nProgram wykorzystuje metode elementow skonczonych oraz rozwiazanie zagadnienia Boussinesqa do wyznaczenia wynikow. Program jest podejsciem do problemu polprzestrzeni od strony  numerycznej. Wykorzystanie Metody Elementow Skonczonych, oraz aproksymacja wspolczynnikow osiadan, sa obarczone bledem metody. Program uzyskuje wyniki w postaci dyskretnej.")


def clicked2():
    messagebox.showinfo("O mnie: inż.Aleksander Mróz",
                        "W czasie gdy pisalem ta wiadomosc, bylem studentem Budownictwa Ladowego i Wodnego na Politechnice Wroclawskiej. Jako specjalizację na studiach magisterskich wybralem Teorie Konstrukcji.")


def clicked3():
    messagebox.showwarning("Nie spodziewałem się hiszpańskiej inkwizicji",
                           "Nikt nie spodziewa się hiszpańskiej inkwizycji")


def clicked4():
    messagebox.showwarning("Poradnik",
                           "Drogi Użytkowniku, pelna wersja instrukcji jest zamieszczona w tresci mojej pracy magisterskiej. Znajdziesz tam pelen opis algorytmu. \n"
                           "\n- Najwazniejsze uwagi:"
                           "\n- Podzial plyty na elementy obliczeniowe nalezy przeprowadzic z glowa. Na boku musza byc min. 4 elementy. \nStosunek bokow elementow nie powinien przekraczac  wartosci 5:1. "
                           "\n- Zanim przeprowadzisz obliczenia, musisz wprowadzić chociaz jedna sile skupiona."
                           "\n- Nie jest wymagane wprowadzanie deformacji gorniczych, aby przeprowadzic obliczenia"
                           "\n\n\n- JESLI cos nie dziala, to nalezy sprawdzic wprowadzone dane. Separator w liczbach w zapisie dziesietnym jest kropka, NIE przecinek"
                           "\n- Albo dziala, ale sie jeszcze liczy. Program lubi w czasie obliczen przejsc w tryb 'brak odpowiedzi'. Obliczenia wykonuja sie w tle.\n"
                           "- Zlozonosc algorytmu jest 2n^3 gdzie n to liczba elementow na ktore zostal podzielony element. Szacunowy czas obliczen dla 400 elementow to okolo 5min"
                           "\n- Przed przeprowadzeniem nalezy poprawnie zdefiniowac poprawnie sciezke do folderu")


class PageObciążenia(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Obciążenia!!!", font=LARGE_FONT)
        # label.grid(column=0, row=1)
        button1 = tk.Button(self, text="Do strony głównej",
                            command=lambda: controller.show_frame(StartPage))
        button1.grid(column=0, row=1)

        #########################################################################

        labe = tk.Label(self, text="Obciążenia\n 'skupione'", font=LARGE_FONT)
        labe.grid(column=1, row=2)
        label1 = tk.Label(self, text="x[m]:", font=LARGE_FONT)
        label1.grid(column=0, row=3)
        Tx = tk.Entry(self, width=10)
        Tx.grid(column=1, row=3)

        label2 = tk.Label(self, text="y[m]:", font=LARGE_FONT)
        label2.grid(column=0, row=4)
        Ty = tk.Entry(self, width=10)
        Ty.grid(column=1, row=4)

        label3 = tk.Label(self, text="F[MN]:", font=LARGE_FONT)
        label3.grid(column=0, row=5)
        Tf = tk.Entry(self, width=10)
        Tf.grid(column=1, row=5)
        label4 = tk.Label(self, text="(Obciążenia zgodne \nz grawitacją => +", font=LARGE_FONT)
        label4.grid(column=1, row=6)

        button2 = tk.Button(self, text="Dodaj Obciążenie",
                            command=lambda: DodajOBC(Tx.get(), Ty.get(), Tf.get()))
        button2.grid(column=1, row=7)
        ###########################################################################################
        przerwa1 = tk.Label(self, text="      ", font=LARGE_FONT)
        przerwa1.grid(column=2, row=1)
        ###########################################################################################
        l0 = tk.Label(self, text="Deformacje \ngórnicze", font=LARGE_FONT)
        l0.grid(column=3, row=2)
        l1 = tk.Label(self, text="R[m]:", font=LARGE_FONT)
        l1.grid(column=2, row=3)
        t1 = tk.Entry(self, width=10)
        t1.grid(column=3, row=3)

        l2x = tk.Label(self, text="Translacja ukladu:", font=LARGE_FONT)
        l2x.grid(column=2, row=4)
        l2 = tk.Label(self, text="DX:", font=LARGE_FONT)
        l2.grid(column=2, row=5)
        t2 = tk.Entry(self, width=10)
        t2.grid(column=3, row=5)

        l3 = tk.Label(self, text="DY:", font=LARGE_FONT)
        l3.grid(column=2, row=6)
        t3 = tk.Entry(self, width=10)
        t3.grid(column=3, row=6)

        przerwa2 = tk.Label(self, text="      ", font=LARGE_FONT)
        przerwa2.grid(column=4, row=0)

        l5 = tk.Label(self, text="kat obrotu [rad]", font=LARGE_FONT)
        l5.grid(column=2, row=7)
        t5 = tk.Entry(self, width=10)
        t5.grid(column=3, row=7)

        przerwa3 = tk.Label(self, text="   ", font=LARGE_FONT)
        przerwa3.grid(column=7, row=0)

        button4 = tk.Button(self, text="Dodaj Front\nGórniczy",
                            command=lambda: DodajFront(t1.get(), t2.get(), t3.get(), t5.get()))
        button4.grid(column=3, row=8)

        Sily = tk.Label(self, text="Tabela przyłożonych\n obciążeń", font=LARGE_FONT)
        Sily.grid(column=1, row=10)

        Gor = tk.Label(self, text="Tabela przyłożonych\n deformacji górniczych", font=LARGE_FONT)
        Gor.grid(column=3, row=10)

        def DodajOBC(*args):
            print(*args)
            X, Y, F = args
            D = [float(X), float(Y), float(F)]
            controller.F.append(D)
            print(controller.F, type(X))
            Var = Sily.cget("text")
            Var += "\n x:,y:,F:   "
            Var += X
            Var += "m , "
            Var += Y
            Var += "m , "
            Var += F
            Var += "MN"
            Sily.configure(text=Var)

        def DodajFront(*args):
            print(*args)
            R, XXX, YYY, alpha = args
            T = [R, XXX, YYY, alpha]
            controller.D.append(T)
            print(controller.D)
            Var = Gor.cget("text")
            Var += "\n X:,Y:,R:,kąt obrotu:   "
            Var += XXX
            Var += "m , "
            Var += YYY
            Var += "m , "
            Var += R
            Var += "m , "
            Var += alpha
            Var += "rad "

            Gor.configure(text=Var)


class HELPOBC(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self,
                         text="Ilustracja Graficzna\n Na czerwowno - deformacje górnicze \n Na zielono - siły skupione",
                         font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Powrót",
                             command=lambda: controller.show_frame(PageObciążenia))
        button1.pack()

        ## image = plt.imread('InstrukcjaOBC.png')
        fig = plt.figure(figsize=(6, 4))
        ## im = plt.imshow(image)  # later use a.set_data(new_data)
        ax = plt.gca()
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


class HELPSTART(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Ilustracja graficzna", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Powrót",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()

        image = plt.imread('MAIN.png')
        fig = plt.figure(figsize=(4, 3))
        im = plt.imshow(image)  # later use a.set_data(new_data)
        ax = plt.gca()
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


class PageMXX(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Wykres Momentów MXX\n wyniki w [MNm]", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Do Strony Głównej",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Pokaż Wykres",
                             command=lambda: self.Pokaz_wykres(controller.MXX, controller.Krotka))
        button2.pack()

    def Pokaz_wykres(self, *args):
        K = args[1]
        a0 = float(K[0])
        b0 = float(K[1])
        m = int(K[2])
        n = int(K[3])
        print("KROTKAwWYK", a0, b0, m, n)
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Make data.
        X = np.arange(a0 / (2 * m), a0, a0 / m)  # Czyli siatkowanie

        Y = np.arange(b0 / (2 * n), b0, b0 / n)  # Siatkowanie

        print(X, Y)

        X, Y = np.meshgrid(X, Y)

        Z = args[0]
        P, O = np.shape(Z)
        min = 9999999999999999999
        max = -999999999999999999
        for i in range(P):
            for j in range(O):
                if (Z[i][j] > max):
                    max = Z[i][j]
                if (Z[i][j] < min):
                    min = Z[i][j]

        print(np.shape(Z), np.shape(X), np.shape(Y))

        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(min - 0.5 * min, max * 1.5)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)

        canvas = FigureCanvasTkAgg(fig, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


class PageMYY(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Wykres momentów MYY\n wyniki w [MNm]", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Do Strony Głównej",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Pokaż Wykres",
                             command=lambda: self.Pokaz_wykres(controller.MYY, controller.Krotka))
        button2.pack()

    def Pokaz_wykres(self, *args):
        K = args[1]
        a0 = float(K[0])
        b0 = float(K[1])
        m = int(K[2])
        n = int(K[3])
        print("KROTKAwWYK", a0, b0, m, n)
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Make data.
        X = np.arange(a0 / (2 * m), a0, a0 / m)  # Czyli siatkowanie

        Y = np.arange(b0 / (2 * n), b0, b0 / n)  # Siatkowanie

        print(X, Y)

        X, Y = np.meshgrid(X, Y)

        Z = args[0]
        P, O = np.shape(Z)
        min = 9999999999999999999
        max = -999999999999999999
        for i in range(P):
            for j in range(O):
                if (Z[i][j] > max):
                    max = Z[i][j]
                if (Z[i][j] < min):
                    min = Z[i][j]

        print(np.shape(Z), np.shape(X), np.shape(Y))

        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(min - 0.5 * min, max * 1.5)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)

        canvas = FigureCanvasTkAgg(fig, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


class PageMXY(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Wykres momentów MXY\n wyniki w [MNm]", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Do Strony Głównej",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Pokaż Wykres",
                             command=lambda: self.Pokaz_wykres(controller.MXY, controller.Krotka))
        button2.pack()

    def Pokaz_wykres(self, *args):
        K = args[1]
        a0 = float(K[0])
        b0 = float(K[1])
        m = int(K[2])
        n = int(K[3])
        print("KROTKAwWYK", a0, b0, m, n)
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Make data.
        X = np.arange(a0 / (2 * m), a0, a0 / m)  # Czyli siatkowanie

        Y = np.arange(b0 / (2 * n), b0, b0 / n)  # Siatkowanie

        print(X, Y)

        X, Y = np.meshgrid(X, Y)

        Z = args[0]
        P, O = np.shape(Z)
        min = 9999999999999999999
        max = -999999999999999999
        for i in range(P):
            for j in range(O):
                if (Z[i][j] > max):
                    max = Z[i][j]
                if (Z[i][j] < min):
                    min = Z[i][j]

        print(np.shape(Z), np.shape(X), np.shape(Y))

        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(min - 0.5 * min, max * 1.5)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)

        canvas = FigureCanvasTkAgg(fig, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


class PageR(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Wykres sił Odporu\n Wyniki w [MN]", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Do Strony Głównej",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Pokaż Wykres",
                             command=lambda: self.Pokaz_wykres(controller.R, controller.Krotka))
        button2.pack()

    def Pokaz_wykres(self, *args):
        K = args[1]
        a0 = float(K[0])
        b0 = float(K[1])
        m = int(K[2])
        n = int(K[3])
        print("KROTKAwWYK", a0, b0, m, n)
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Make data.
        X = np.arange(a0 / (2 * m), a0, a0 / m)  # Czyli siatkowanie

        Y = np.arange(b0 / (2 * n), b0, b0 / n)  # Siatkowanie

        print(X, Y)

        X, Y = np.meshgrid(X, Y)

        Z = args[0]
        P, O = np.shape(Z)
        min = 9999999999999999999
        max = -999999999999999999
        for i in range(P):
            for j in range(O):
                if (Z[i][j] > max):
                    max = Z[i][j]
                if (Z[i][j] < min):
                    min = Z[i][j]

        print(np.shape(Z), np.shape(X), np.shape(Y))

        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(min - 0.5 * min, max * 1.5)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)

        canvas = FigureCanvasTkAgg(fig, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


class PageOsiadania(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Wykres Osiadań\n wyniki w [m]", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Do Strony Głównej",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Pokaż Wykres",
                             command=lambda: self.Pokaz_wykres(controller.W, controller.Krotka))
        button2.pack()

    def Pokaz_wykres(self, *args):
        K = args[1]
        a0 = float(K[0])
        b0 = float(K[1])
        m = int(K[2])
        n = int(K[3])
        print("KROTKAwWYK", a0, b0, m, n)
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Make data.
        X = np.arange(a0 / (2 * m), a0, a0 / m)  # Czyli siatkowanie

        Y = np.arange(b0 / (2 * n), b0, b0 / n)  # Siatkowanie

        print(X, Y)

        X, Y = np.meshgrid(X, Y)

        Z = args[0]
        P, O = np.shape(Z)
        min = 9999999999999999999
        max = -999999999999999999
        for i in range(P):
            for j in range(O):
                if (Z[i][j] > max):
                    max = Z[i][j]
                if (Z[i][j] < min):
                    min = Z[i][j]

        print(np.shape(Z), np.shape(X), np.shape(Y))

        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(min - 0.5 * min, max * 1.5)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)

        canvas = FigureCanvasTkAgg(fig, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


def animate(i):
    pullData = open('sampleText.txt', 'r').read()

    dataArray = pullData.split('\n')

    xar = []
    yar = []
    for eachLine in dataArray:
        if len(eachLine) > 1:
            x, y = eachLine.split(',')
            xar.append(int(x))
            yar.append(int(y))
    a.clear()
    a.plot(xar, yar)


app = SeaofBTCapp()

app.mainloop()
