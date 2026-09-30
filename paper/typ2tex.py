"""Convert paper.typ (the subset of Typst used there) to LaTeX. Heuristic but complete for this document."""
import re, sys

src = open('paper.typ').read()

# ---------- math conversion
SYM = {
    'RR': r'\mathbb{R}', 'ZZ': r'\mathbb{Z}', 'without': r'\setminus', 'arrow': r'\to', 'arrow.bar': r'\mapsto',
    'subset.eq': r'\subseteq', 'subset': r'\subset', 'union.big': r'\bigcup', 'union': r'\cup', 'inter': r'\cap',
    'chevron.l': r'\langle', 'chevron.r': r'\rangle', 'dot': r'\cdot', 'times.r': r'\rtimes', 'times': r'\times',
    'plus.minus': r'\pm', 'infinity': r'\infty', 'emptyset': r'\emptyset', 'partial': r'\partial',
    'epsilon': r'\varepsilon', 'gamma': r'\gamma', 'Gamma': r'\Gamma', 'Sigma': r'\Sigma', 'pi': r'\pi',
    'sigma': r'\sigma', 'lambda': r'\lambda', 'in': r'\in', 'dots': r'\dots', 'sup': r'\sup', 'sum': r'\sum',
    'ker': r'\ker', 'id': r'\mathrm{id}', 'quad': r'\quad', 'approx': r'\approx', 'square': r'\square',
    'nothing': r'\emptyset', 'infty': r'\infty',
}

def conv_math(m):
    s = m
    s = s.replace('\\/', '/')
    # quoted text
    s = re.sub(r'"([^"]*)"', lambda k: '\\mathrm{' + k.group(1).replace(' ', '\\ ') + '}', s)
    # cal(X)
    s = re.sub(r'cal\(([A-Za-z])\)', r'\\mathcal{\1}', s)
    # multi-letter symbols (longest first)
    for k in sorted(SYM, key=len, reverse=True):
        s = re.sub(r'(?<![A-Za-z\\])' + re.escape(k) + r'(?![A-Za-z.])', lambda _m, v=SYM[k]: v, s)
    # sub/superscripts with parentheses: _(...) -> _{...}
    def grp(s, ch):
        out = ''; i = 0
        while i < len(s):
            if s[i] == ch and i + 1 < len(s) and s[i+1] == '(':
                depth = 0; j = i + 1
                while j < len(s):
                    if s[j] == '(': depth += 1
                    elif s[j] == ')':
                        depth -= 1
                        if depth == 0: break
                    j += 1
                inner = s[i+2:j]
                out += ch + '{' + grp(inner, ch) + '}'; i = j + 1
            else:
                out += s[i]; i += 1
        return out
    s = grp(grp(s, '_'), '^')
    # multi-char sub/superscripts without parens e.g. ^2 fine; x_ab -> only first char in typst; leave.
    s = s.replace('>=', r'\ge ').replace('<=', r'\le ').replace('!=', r'\ne ')
    s = s.replace('|', '|')
    s = re.sub(r'\bN 2\^', r'N\\,2^', s)
    return s

def conv_inline_math(text):
    def repl(m):
        inner = m.group(1)
        if inner.startswith(' ') and inner.endswith(' '):
            return '\\[' + conv_math(inner.strip()) + '\\]'
        return '$' + conv_math(inner) + '$'
    return re.sub(r'\$(.+?)\$', repl, text, flags=re.S)

# ---------- text conversion helpers
def esc(t):
    return t

def conv_text(t):
    t = conv_inline_math(t)
    # protect math from the markup passes below
    stash = []
    def _st(m):
        stash.append(m.group(0)); return '\x00%d\x00' % (len(stash)-1)
    t = re.sub(r'\\\[.*?\\\]|\$[^$]*\$', _st, t, flags=re.S)
    t = _conv_text_rest(t)
    return re.sub('\x00(\\d+)\x00', lambda m: stash[int(m.group(1))], t)

def _conv_text_rest(t):
    # citations (only keys not followed by ':')
    t = re.sub(r'((?:@[a-z0-9]+(?![a-z0-9:])\s*)+)', lambda m: r'\cite{' + ','.join(re.findall(r'@([a-z0-9]+)', m.group(1))) + '}' + (' ' if m.group(1).endswith(' ') else ''), t)
    # refs
    t = re.sub(r'@sec:([a-zA-Z0-9]+)', r'Section~\\ref{sec:\1}', t)
    t = re.sub(r'@fig:([a-zA-Z0-9]+)', r'Figure~\\ref{fig:\1}', t)
    t = re.sub(r'@tab:([a-zA-Z0-9]+)', r'Table~\\ref{tab:\1}', t)
    # bold / italic (outside math)
    t = re.sub(r'\*([^*\n]+)\*', r'\\textbf{\1}', t)
    t = re.sub(r'(?<![A-Za-z0-9\\{])_([^_\n]+?)_(?![A-Za-z0-9])', r'\\emph{\1}', t)
    # raw code
    t = re.sub(r'`([^`]+)`', lambda m: r'\texttt{' + m.group(1).replace('_', r'\_') + '}', t)
    # links
    t = re.sub(r'#link\("([^"]+)"\)\[([^\]]+)\]', r'\\href{\1}{\2}', t)
    t = t.replace('—', '---').replace('–', '--').replace('#h(1fr)', '')
    t = t.replace('&', r'\&').replace('#', r'\#') if False else t
    return t

lines = src.split('\n')
out = []
i = 0
# skip preamble until title
while i < len(lines) and not lines[i].startswith('#align(center)['): i += 1
i += 1
while not lines[i].startswith(']'): i += 1
i += 1

def find_block_end(lines, i):
    """given line index of a line containing '[' opening a content block that ends at a line '...]' return index."""
    depth = 0; j = i
    while j < len(lines):
        depth += lines[j].count('[') - lines[j].count(']')
        if depth <= 0: return j
        j += 1
    return j

body = '\n'.join(lines[i:])

# abstract block
m = re.search(r'#block\(inset: \(x: 1\.5em\)\)\[\n\*Abstract\.\* (.*?)\n\n_Author contribution\._ (.*?)\n\]', body, re.S)
abstract = conv_text(m.group(1)); contribution = conv_text(m.group(2))
body = body[m.end():]

# theorem environments: #thm("Kind", [name])[ ... ]  possibly multi-line; #pf[ ... ]

def match_bracket(body, j):
    """j = index just after the opening '['; return index just after the matching ']' (skips $...$)."""
    depth = 1
    while depth > 0:
        c = body[j]
        if c == '$':
            j = body.index('$', j+1) + 1; continue
        if c == '[': depth += 1
        elif c == ']': depth -= 1
        j += 1
    return j

def replace_thm(body):
    res = ''; pos = 0
    for m in re.finditer(r'#thm\("([^"]+)", (\[[^\]]*\]|none)\)\[', body):
        if m.start() < pos: continue
        res += body[pos:m.start()]
        # find matching bracket for the body
        j = match_bracket(body, m.end())
        content = body[m.end():j-1]
        kind = m.group(1); name = m.group(2)
        name = None if name == 'none' else name[1:-1]
        env = {'Definition': 'definition', 'Lemma': 'lemma', 'Theorem': 'theorem', 'Conjecture': 'conjecture',
               'Construction sketch': 'construction'}.get(kind)
        label = ''
        if kind.startswith('Proposition'): env = 'proposition'; label = r'\label{prop:cert}'
        elif kind.startswith('Problem'): env = 'problem'; label = r'\label{prob:N4}'
        elif kind.startswith('Theorem 11'): env = 'theorem'; label = r'\label{thm:cond}'
        elif kind == 'Theorem A': env = 'theoremA'
        elif kind == 'Theorem B': env = 'theoremB'
        opt = ('[' + name + ']') if name else ''
        res += '\\begin{' + env + '}' + opt + label + '\n' + content.strip() + '\n\\end{' + env + '}\n'
        pos = j
    return res + body[pos:]
body = replace_thm(body)
body = re.sub(r'#pf\[', r'\\begin{proof}\n', body)
# close proofs: the #pf[...] closing bracket is the first ']' at line end after; handle by matching
def close_proofs(body):
    res = ''; pos = 0
    while True:
        k = body.find('\\begin{proof}\n', pos)
        if k < 0: break
        j = match_bracket(body, k + len('\\begin{proof}\n'))
        res += body[pos:j-1] + '\n\\end{proof}\n'; pos = j
    return res + body[pos:]
body = close_proofs(body)

# figures with images
def fig_repl(m):
    path, width, cap, lab = m.group(1), m.group(2), m.group(3), m.group(4)
    return ('\\begin{figure}[t]\\centering\\includegraphics[width=' + width + '\\linewidth]{' + path + '}\n'
            '\\caption{' + cap + '}\\label{' + lab + '}\\end{figure}')
body = re.sub(r'#figure\(image\("([^"]+)", width: (\d+)%\), caption: \[(.*?)\]\) <([a-z:A-Z0-9]+)>',
              lambda m: fig_repl(type('M', (), {'group': lambda self, k: [None, m.group(1), str(int(m.group(2))/100), m.group(3), m.group(4)][k]})()), body, flags=re.S)

# tables
def table_repl(m):
    head, rows, cap, lab = m.group(1), m.group(2), m.group(3), m.group(4)
    ncol = len(re.findall(r'columns: \(([^)]*)\)', head)[0].split(',')) if 'columns: (' in head else int(re.findall(r'columns: (\d+)', head)[0])
    cells = re.findall(r'\[((?:[^\[\]]|\[[^\]]*\])*)\]', rows)
    trs = []
    for r in range(0, len(cells), ncol):
        trs.append(' & '.join(cells[r:r+ncol]) + r' \\')
    colspec = 'l' * ncol if ncol != 6 else 'cccccc'
    if ncol == 6: colspec = 'cp{2.6cm}ccc p{4.3cm}'
    if ncol == 4: colspec = 'p{5.2cm}rp{3.2cm}p{4cm}'
    if ncol == 2 and 'quantity' in rows: colspec = 'p{4.5cm}p{8cm}'
    tex = '\\begin{table}[t]\\centering\\small\n\\begin{tabular}{' + colspec + '}\\toprule\n' + trs[0] + ' \\midrule\n' + '\n'.join(trs[1:]) + '\n\\bottomrule\\end{tabular}\n\\caption{' + cap + '}\\label{' + lab + '}\\end{table}'
    return tex
body = re.sub(r'#figure\(table\(([^\n]*)\n(.*?)\n\), caption: \[(.*?)\]\) <([a-z:A-Z0-9]+)>', table_repl, body, flags=re.S)

# headings
body = re.sub(r'^= (.*?)(?: <sec:([a-zA-Z0-9]+)>)?$', lambda m: '\\section{' + m.group(1) + '}' + ('\\label{sec:' + m.group(2) + '}' if m.group(2) else ''), body, flags=re.M)
body = re.sub(r'#heading\(numbering: none\)\[(.*?)\]', r'\\section*{\1}', body)
# lists
body = re.sub(r'((?:^\d+\. .*\n(?:\n(?=\d+\. ))?)+)', lambda m: '\\begin{enumerate}\n' + re.sub(r'^\d+\. ', r'\\item ', m.group(1), flags=re.M) + '\\end{enumerate}\n', body, flags=re.M)
body = re.sub(r'((?:^- .*\n(?:\n(?=- ))?)+)', lambda m: '\\begin{itemize}\n' + re.sub(r'^- ', r'\\item ', m.group(1), flags=re.M) + '\\end{itemize}\n', body, flags=re.M)
# appendix A poses grid
_poses = '\\begin{multicols}{3}\\scriptsize\\begin{verbatim}\n' + open('poses44.txt').read() + '\n\\end{verbatim}\\end{multicols}'
body = re.sub(r'#text\(7\.5pt\)\[#grid\(.*?\)\]', lambda m: _poses, body, flags=re.S)
# bibliography
body = body.replace('#bibliography("refs.yml", title: "References", style: "ieee")', '\\bibliographystyle{plain}\n\\bibliography{refs}')
body = conv_text(body)
# textbf inside theorem names etc fine; fix Theorem A/B references
body = body.replace('Proposition 7.1', 'Proposition~\\ref{prop:cert}').replace('Problem 9.1', 'Problem~\\ref{prob:N4}').replace('Theorem 11.1', 'Theorem~\\ref{thm:cond}')
body = body.replace('\\begin{proposition}[finite certificate]', '\\begin{proposition}[finite certificate]')
body = body.replace('<sec:contrib>', '').replace('#v(0.3em)', '')
body = body.replace('%', r'\%')
# citation inside \cite after conversion may include labels like tsiokos; ok.

preamble = r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{multicol}
\usepackage[margin=2.4cm]{geometry}

\usepackage{lmodern}
\usepackage{microtype}
\usepackage[hidelinks]{hyperref}
\theoremstyle{plain}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{conjecture}[theorem]{Conjecture}
\newtheorem{problem}[theorem]{Problem}
\newtheorem*{theoremA}{Theorem A}
\newtheorem*{theoremB}{Theorem B}
\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{construction}[theorem]{Construction sketch}
\title{Towards Strongly Aperiodic Monotiles\\ in Higher Dimensions}
\author{Dmitry Kamenetsky\\ \texttt{dkamenetsky@gmail.com}}
\date{September 30, 2026}
\begin{document}
\maketitle
\begin{abstract}
''' + abstract + r'''
\end{abstract}
\noindent\emph{Author contribution.} ''' + contribution + '\n\n'

open('paper.tex', 'w').write(preamble + body + '\n\\end{document}\n')
print('written')
