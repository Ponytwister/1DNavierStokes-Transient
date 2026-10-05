# Application icon

`navier-source.png` is the full-resolution artwork selected as Option C, generated
with the built-in image generation tool from the original microchannel icon.
It depicts a pale-blue chip with three navy inlet channels merging into one
outlet and a coral-red detection zone. The transparent source is preserved for
regeneration; the previous isometric SVG is available in Git history.

`navier.ico` contains transparent 16, 20, 24, 32, 40, 48, 64, 128, and 256 pixel
versions of this shared application identity. `navier.png` is a 256 pixel preview.
Windows resources embed the icon in both Navier.exe and NavierGui.exe; the Qt
resource also supplies the desktop application's window/taskbar icon.
Test utilities retain their default icons.

The committed assets require no image tools at build time. To regenerate,
downsample `navier-source.png` with Lanczos to the listed sizes and encode them
as a multi-image Windows ICO (RGBA). Keep the PNG preview at 256 pixels and
preserve the source alpha channel. No numerical behavior changes are involved.

Generation prompt: Create one alternative Windows application icon for a
microchannel flow simulation program, using the existing icon as subject and
color reference. Very minimal modern app icon, viewed mostly from above: a
compact pale-blue rounded square microfluidic chip with three thick navy inlet
paths merging into a single outlet, a bright coral-red detection zone crossing
the outlet. Flat geometric shapes, strong negative space, no floating plate.
Bold instantly readable scientific instrument symbol at 16-32 pixels. Square
composition, single centered icon filling about 85% of canvas, genuinely
transparent background with alpha. No text, letters, labels, arrows, watermark,
border, or comparison sheet. Variant C - Minimal channel mark.
