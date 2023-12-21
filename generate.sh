python3 md2nb.py index.md
jupyter nbconvert index.ipynb --to slides --reveal-prefix reveal.js
--SlidesExporter.reveal_theme=black 
--SlidesExporter.reveal_scroll=True 
--SlidesExporter.reveal_transition=none
mv index.slides.html index.html