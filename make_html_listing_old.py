#!/bin/env python3

import os
import sys

if __name__ == '__main__':
    plot_dir = sys.argv[1]
    #html_dir = f'/publicweb/m/mreid/iDM_AN_Plots/{plot_dir}'
    html_dir = f'{plot_dir}'
    if os.path.isdir(html_dir) == False:
        os.mkdir(html_dir)

    html = open(os.path.join(html_dir, 'index.html'), 'wt')
    html.write('<html><head></head><body><pre>')
    file_list = os.listdir(html_dir)
    truncated_file_list = [os.path.splitext(os.path.basename(f))[0] for f in file_list]
    html.write('<a href="..">.. (parent directory)</a><br>')
    for i, f in enumerate(file_list):
        print(i, f)
        if f != 'index.html':
            html.write(f'<a href="#{truncated_file_list[i]}">{i}</a>')
            html.write(f'   ')
            html.write(f'<a href="{f}">{f}</a><br>')

    for i, f in enumerate(file_list):
        if f != 'index.html':
            html.write(f'<h4 id="{truncated_file_list[i]}"><a href="#{truncated_file_list[i]}">{truncated_file_list[i]}</a></h4><br>')
            html.write(f'<img src="{f}" style="max-width: 600px"><br><br>\n')
    html.write('</pre></body></html>\n')
