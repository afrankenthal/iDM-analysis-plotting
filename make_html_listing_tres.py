#!/bin/env python3
import sys
from pathlib import Path

if __name__ == '__main__':
    plot_dir = Path(sys.argv[1])
    html_dir = Path(f'{plot_dir}')
    #html_dir = Path(f'/publicweb/a/as2872/iDMPlots/{plot_dir}')
    if (not html_dir.exists()):
        html_dir.mkdir()

    html_file = Path(html_dir/'index.html')
    with html_file.open('wt') as index:
        index.write('<html><head></head><body><pre>')
        index.write('<a href="..">.. (parent directory)</a><br>')
        for i, f in enumerate(sorted(html_dir.glob('*'))):
            if f.name != 'index.html':
                print(i, f.name)
                index.write(f'<a href="#{f.name}">{i}</a>')
                index.write(f'   ')
                index.write(f'<a href="{f.name}">{f.name}</a><br>')

        for i, f in enumerate(sorted(html_dir.glob('*'))):
            if f.name != 'index.html':
                index.write(f'<h4 id="{f.name}"><a href="#{f.name}">{f.name}</a></h4><br>')
                if f.suffix == '.png':
                    index.write(f'<a href="{f.name}"><img src="{f.name}" style="max-width: 600px"></a><br><br>\n')
                else:
                    index.write(f'<a href="{f.name}">{f.name}</a><br><br>\n')
        index.write('</pre></body></html>\n')
