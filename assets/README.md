# Application icon

`navier.svg` is the editable vector artwork, inspired by the user-supplied
Microchannel.png: a pale-blue microfluidic chip, converging inlet channels, and
a red fluorescence detector outline. Labels and leader arrows are omitted for
legibility. The original reference image is not needed to build the application.

`navier.ico` contains transparent 16, 20, 24, 32, 40, 48, 64, 128, and 256 pixel
versions of this shared application identity. `navier.png` is a 256 pixel preview.
Windows resources embed the icon in both Navier.exe and NavierGui.exe; the Qt
resource also supplies the desktop application's window/taskbar icon.
Test utilities retain their default icons.

The committed assets require no image tools at build time. To regenerate, render
the SVG at 1024 pixels with an SVG renderer, then downsample with Lanczos to the
listed sizes and encode them as a multi-image Windows ICO (RGBA). Keep the PNG
preview at 256 pixels. No runtime or numerical behavior changes are involved.
