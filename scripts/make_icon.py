#!/usr/bin/env python3
"""Generate macOS .icns icon from a source SVG or PNG."""

import os
import subprocess
import sys


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.normpath(os.path.join(script_dir, ".."))
    iconset_dir = os.path.join(repo_root, "resources", "AppIcon.iconset")
    os.makedirs(iconset_dir, exist_ok=True)

    try:
        from PIL import Image, ImageDraw, ImageFont

        sizes = [16, 32, 64, 128, 256, 512, 1024]

        for size in sizes:
            img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
            draw = ImageDraw.Draw(img)

            margin = int(size * 0.1)
            radius = int(size * 0.18)
            draw.rounded_rectangle(
                [margin, margin, size - margin, size - margin],
                radius=radius,
                fill=(20, 60, 90, 255),
            )

            font_size = int(size * 0.55)
            try:
                font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Times New Roman.ttf", font_size)
            except (IOError, OSError):
                try:
                    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf", font_size)
                except (IOError, OSError):
                    font = ImageFont.load_default()

            text = "ψ"
            bbox = draw.textbbox((0, 0), text, font=font)
            text_width = bbox[2] - bbox[0]
            text_height = bbox[3] - bbox[1]
            x = (size - text_width) // 2
            y = (size - text_height) // 2 - int(size * 0.05)
            draw.text((x, y), text, fill=(255, 255, 255, 255), font=font)

            img_resized = img.resize((size, size), Image.LANCZOS)
            img_resized.save(os.path.join(iconset_dir, f"icon_{size}x{size}.png"))

            if size <= 512:
                img_2x = img.resize((size * 2, size * 2), Image.LANCZOS)
                img_2x.save(os.path.join(iconset_dir, f"icon_{size}x{size}@2x.png"))

        print(f"Generated iconset at {iconset_dir}")

    except ImportError:
        print("Pillow not available, creating placeholder icons")
        import struct
        import zlib

        def create_png(width, height, color):
            def chunk(chunk_type, data):
                c = chunk_type + data
                return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)

            header = b"\x89PNG\r\n\x1a\n"
            ihdr = chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0))

            raw_data = b""
            for _ in range(height):
                raw_data += b"\x00"
                for _ in range(width):
                    raw_data += bytes(color)

            idat = chunk(b"IDAT", zlib.compress(raw_data))
            iend = chunk(b"IEND", b"")
            return header + ihdr + idat + iend

        for size in [16, 32, 128, 256, 512]:
            png_data = create_png(size, size, [20, 60, 90, 255])
            with open(os.path.join(iconset_dir, f"icon_{size}x{size}.png"), "wb") as f:
                f.write(png_data)
            png_data_2x = create_png(size * 2, size * 2, [20, 60, 90, 255])
            with open(os.path.join(iconset_dir, f"icon_{size}x{size}@2x.png"), "wb") as f:
                f.write(png_data_2x)

    icns_path = os.path.join(repo_root, "resources", "AppIcon.icns")
    if sys.platform == "darwin":
        subprocess.run(["iconutil", "-c", "icns", iconset_dir, "-o", icns_path], check=True)
        print(f"Generated {icns_path}")
    else:
        print("Not on macOS - skipping .icns generation")


if __name__ == "__main__":
    main()
