"""
Drag-and-drop layout editor for TRENCH faceplate.
Run: python tools/layout_editor.py
Drag the colored rectangles to position them on the faceplate.
Close the window to print the final coordinates.
"""
import tkinter as tk
from PIL import Image, ImageTk

FACEPLATE = "ui_faceplate.png"

# Current layout — edit these starting positions
elements = {
    "MORPH_STRIP":   {"x": 30,  "y": 174, "w": 137, "h": 37, "color": "cyan"},
    "MORPH_READOUT": {"x": 176, "y": 183, "w": 57,  "h": 19, "color": "lime"},
    "Q_STRIP":       {"x": 32,  "y": 237, "w": 132, "h": 31, "color": "cyan"},
    "Q_READOUT":     {"x": 177, "y": 245, "w": 53,  "h": 17, "color": "lime"},
    "TYPE_BAR":      {"x": 70,  "y": 30,  "w": 200, "h": 14, "color": "orange"},
}

class LayoutEditor:
    def __init__(self, root):
        self.root = root
        root.title("TRENCH Layout Editor — drag to position, close to print")

        img = Image.open(FACEPLATE)
        self.tk_img = ImageTk.PhotoImage(img)

        self.canvas = tk.Canvas(root, width=img.width, height=img.height)
        self.canvas.pack()
        self.canvas.create_image(0, 0, anchor=tk.NW, image=self.tk_img)

        self.rects = {}
        self.labels = {}
        self.drag_data = {"item": None, "ox": 0, "oy": 0}

        for name, el in elements.items():
            # Filled rectangle — easy to grab
            r = self.canvas.create_rectangle(
                el["x"], el["y"], el["x"] + el["w"], el["y"] + el["h"],
                outline=el["color"], fill=el["color"], width=2,
                stipple="gray50"
            )
            lbl = self.canvas.create_text(
                el["x"] + el["w"] // 2, el["y"] + el["h"] // 2,
                anchor=tk.CENTER,
                text=name, fill="black", font=("Consolas", 7, "bold")
            )
            self.rects[r] = name
            self.labels[name] = lbl

            self.canvas.tag_bind(r, "<ButtonPress-1>", self.on_press)
            self.canvas.tag_bind(r, "<B1-Motion>", self.on_drag)
            self.canvas.tag_bind(lbl, "<ButtonPress-1>", lambda e, rid=r: self.on_press_redirect(e, rid))
            self.canvas.tag_bind(lbl, "<B1-Motion>", self.on_drag)

        # Coordinate display
        self.coord_var = tk.StringVar(value="Drag a rectangle...")
        tk.Label(root, textvariable=self.coord_var, font=("Consolas", 9),
                 bg="black", fg="white", anchor=tk.W).pack(fill=tk.X)

        root.protocol("WM_DELETE_WINDOW", self.on_close)

    def on_press_redirect(self, event, rect_id):
        self.drag_data["item"] = rect_id
        coords = self.canvas.coords(rect_id)
        self.drag_data["ox"] = event.x - coords[0]
        self.drag_data["oy"] = event.y - coords[1]

    def on_press(self, event):
        item = self.canvas.find_closest(event.x, event.y)[0]
        if item in self.rects:
            self.drag_data["item"] = item
            coords = self.canvas.coords(item)
            self.drag_data["ox"] = event.x - coords[0]
            self.drag_data["oy"] = event.y - coords[1]

    def on_drag(self, event):
        item = self.drag_data["item"]
        if item is None:
            return
        name = self.rects[item]
        el = elements[name]
        nx = event.x - self.drag_data["ox"]
        ny = event.y - self.drag_data["oy"]
        # Update element position
        el["x"] = int(nx)
        el["y"] = int(ny)
        self.canvas.coords(item, nx, ny, nx + el["w"], ny + el["h"])
        # Move label to center of rect
        self.canvas.coords(self.labels[name], nx + el["w"] // 2, ny + el["h"] // 2)
        # Update status
        self.coord_var.set(f"{name}: x={el['x']}, y={el['y']}")

    def on_close(self):
        print("\n// ---- FINAL COORDINATES ---- paste into editor.rs lines 37-63 ----")
        print()
        m = elements["MORPH_STRIP"]
        print(f"const MORPH_STRIP_X: f32 = {m['x']}.0;")
        print(f"const MORPH_STRIP_Y: f32 = {m['y']}.0;")
        print(f"const STRIP_W: f32 = {m['w']}.0;")
        print(f"const STRIP_H: f32 = {m['h']}.0;")
        print()
        q = elements["Q_STRIP"]
        print(f"const Q_STRIP_X: f32 = {q['x']}.0;")
        print(f"const Q_STRIP_Y: f32 = {q['y']}.0;")
        print(f"const Q_STRIP_W: f32 = {q['w']}.0;")
        print(f"const Q_STRIP_H: f32 = {q['h']}.0;")
        print()
        mr = elements["MORPH_READOUT"]
        print(f"const MORPH_READOUT_X: f32 = {mr['x']}.0;")
        print(f"const MORPH_READOUT_Y: f32 = {mr['y']}.0;")
        print(f"const MORPH_READOUT_W: f32 = {mr['w']}.0;")
        print(f"const MORPH_READOUT_H: f32 = {mr['h']}.0;")
        print()
        qr = elements["Q_READOUT"]
        print(f"const Q_READOUT_X: f32 = {qr['x']}.0;")
        print(f"const Q_READOUT_Y: f32 = {qr['y']}.0;")
        print(f"const Q_READOUT_W: f32 = {qr['w']}.0;")
        print(f"const Q_READOUT_H: f32 = {qr['h']}.0;")
        print()
        t = elements["TYPE_BAR"]
        print(f"const TYPE_BAR_X: f32 = {t['x']}.0;")
        print(f"const TYPE_BAR_Y: f32 = {t['y']}.0;")
        print(f"const TYPE_BAR_W: f32 = {t['w']}.0;")
        print(f"const TYPE_BAR_H: f32 = {t['h']}.0;")
        print()
        print("// ---- END ----")
        self.root.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    LayoutEditor(root)
    root.mainloop()
