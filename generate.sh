python3 md2nb.py index.md
jupyter nbconvert index.ipynb --to slides \
--SlidesExporter.reveal_scroll=True
mv index.slides.html index.html