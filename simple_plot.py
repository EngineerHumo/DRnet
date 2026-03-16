import math
import struct
import zlib


class Canvas:
    def __init__(self, width=900, height=650, bg=(255, 255, 255)):
        self.width = int(width)
        self.height = int(height)
        self.px = [list(bg) for _ in range(self.width * self.height)]

    def _idx(self, x, y):
        return y * self.width + x

    def set_pixel(self, x, y, color):
        if 0 <= x < self.width and 0 <= y < self.height:
            self.px[self._idx(int(x), int(y))] = [int(color[0]), int(color[1]), int(color[2])]

    def fill_rect(self, x0, y0, x1, y1, color):
        xa, xb = sorted((int(x0), int(x1)))
        ya, yb = sorted((int(y0), int(y1)))
        xa = max(0, xa)
        ya = max(0, ya)
        xb = min(self.width - 1, xb)
        yb = min(self.height - 1, yb)
        for y in range(ya, yb + 1):
            base = y * self.width
            for x in range(xa, xb + 1):
                self.px[base + x] = [color[0], color[1], color[2]]

    def line(self, x0, y0, x1, y1, color=(0, 0, 0), width=1):
        x0, y0, x1, y1 = int(x0), int(y0), int(x1), int(y1)
        dx = abs(x1 - x0)
        sx = 1 if x0 < x1 else -1
        dy = -abs(y1 - y0)
        sy = 1 if y0 < y1 else -1
        err = dx + dy
        while True:
            for ox in range(-width // 2, width // 2 + 1):
                for oy in range(-width // 2, width // 2 + 1):
                    self.set_pixel(x0 + ox, y0 + oy, color)
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                x0 += sx
            if e2 <= dx:
                err += dx
                y0 += sy

    def circle(self, cx, cy, r=3, color=(0, 0, 0), fill=True):
        cx, cy, r = int(cx), int(cy), int(r)
        for y in range(cy - r, cy + r + 1):
            for x in range(cx - r, cx + r + 1):
                d = (x - cx) * (x - cx) + (y - cy) * (y - cy)
                if (fill and d <= r * r) or (not fill and abs(d - r * r) < r):
                    self.set_pixel(x, y, color)

    def to_png_bytes(self):
        raw = bytearray()
        for y in range(self.height):
            raw.append(0)
            start = y * self.width
            for x in range(self.width):
                raw.extend(bytes(self.px[start + x]))

        def chunk(tag, data):
            body = tag + data
            return struct.pack('!I', len(data)) + body + struct.pack('!I', zlib.crc32(body) & 0xFFFFFFFF)

        out = bytearray(b'\x89PNG\r\n\x1a\n')
        ihdr = struct.pack('!IIBBBBB', self.width, self.height, 8, 2, 0, 0, 0)
        out.extend(chunk(b'IHDR', ihdr))
        out.extend(chunk(b'IDAT', zlib.compress(bytes(raw), level=9)))
        out.extend(chunk(b'IEND', b''))
        return bytes(out)

    def save_png(self, path):
        with open(path, 'wb') as f:
            f.write(self.to_png_bytes())

    def save_pdf(self, path):
        rgb = bytearray()
        for p in self.px:
            rgb.extend(bytes(p))
        comp = zlib.compress(bytes(rgb), 9)
        objs = []
        objs.append(b"1 0 obj<< /Type /Catalog /Pages 2 0 R >>endobj\n")
        objs.append(b"2 0 obj<< /Type /Pages /Kids [3 0 R] /Count 1 >>endobj\n")
        objs.append(f"3 0 obj<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self.width} {self.height}] /Contents 4 0 R /Resources << /XObject << /Im0 5 0 R >> >> >>endobj\n".encode())
        content = f"q\n{self.width} 0 0 {self.height} 0 0 cm\n/Im0 Do\nQ\n".encode()
        objs.append(b"4 0 obj<< /Length " + str(len(content)).encode() + b" >>stream\n" + content + b"endstream\nendobj\n")
        img_hdr = f"5 0 obj<< /Type /XObject /Subtype /Image /Width {self.width} /Height {self.height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length {len(comp)} >>stream\n".encode()
        objs.append(img_hdr + comp + b"\nendstream\nendobj\n")

        out = bytearray(b"%PDF-1.4\n")
        offsets = [0]
        for o in objs:
            offsets.append(len(out))
            out.extend(o)
        xref = len(out)
        out.extend(f"xref\n0 {len(objs)+1}\n".encode())
        out.extend(b"0000000000 65535 f \n")
        for off in offsets[1:]:
            out.extend(f"{off:010d} 00000 n \n".encode())
        out.extend(f"trailer<< /Size {len(objs)+1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode())
        with open(path, 'wb') as f:
            f.write(bytes(out))


def scale(v, lo, hi, out_lo, out_hi):
    if hi - lo == 0:
        return (out_lo + out_hi) / 2
    return out_lo + (v - lo) * (out_hi - out_lo) / (hi - lo)


def quantile(vals, q):
    arr = sorted(vals)
    if not arr:
        return 0.0
    k = (len(arr) - 1) * q
    f = int(math.floor(k))
    c = int(math.ceil(k))
    if f == c:
        return arr[f]
    return arr[f] * (c - k) + arr[c] * (k - f)
