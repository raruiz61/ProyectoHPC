'''
Parses an mae file into a list/dict like datastructure.

Example:

    >>> reader = M2IOReader('hypothesis.xvol')
    >>> block = reader.find_block('f_m_table')
    >>> print 'radii', block.m_row.r_m_phase_radius

'''

import re

re_m = re.compile(r'(\w+)\[(\d+)\]$')
re_space = re.compile(r'\s')

# maps the type prefix of column names to cast callables
typemap = {
    'i': lambda i: 0 if i == '<>' else int(i),
    'r': lambda i: 0.0 if i == '<>' else float(i),
    's': str,
}
typemap['b'] = typemap['i']

def m2io_repr(value, key=None):
    '''
    M2IO string of basic types or types with an 'm2io_repr' method.
    '''
    try:
        return value.m2io_repr(key)
    except AttributeError:
        pass
    r = str(value)
    if re_space.search(r):
        r = '"' + r.replace('"', r'\"') + '"'
    return r

class M2IO(list):
    '''
    List-like datastructure of M2IO blocks.
    '''
    def m2io_repr(self, name=''):
        return '\n'.join(v.m2io_repr() for v in self)

class Table(list):
    '''
    M2IO Table.

    instance.keys -> list of column names
    instance[0] -> first row
    instance.abc -> column with name "abc"
    '''
    def __init__(self, keys):
        self._types = [typemap[k[0]] for k in keys]
        self.keys = keys

    def add_row(self, iterable):
        '''
        Append a data row to the end of the table
        '''
        self.append([t(c) for (t, c) in zip(self._types, iterable)])

    def __repr__(self):
        return 'Table(<%d>)' % len(self)

    def __getattr__(self, key):
        i = self.keys.index(key)
        return [r[i] for r in self]

    def m2io_repr(self, name=''):
        r = '%s[%d] {\n' % (name, len(self))
        r += ' '.join(self.keys)
        r += '\n:::\n'
        for i, line in enumerate(self, 1):
            r += str(i)
            for value in line:
                r += ' ' + m2io_repr(value)
            r += '\n'
        r += ':::\n}\n'
        return r

class Block(dict):
    '''
    M2IO top level block
    '''
    def __init__(self, name, items):
        dict.__init__(self, items)
        self.name = name

    def __repr__(self):
        r = dict.__repr__(self)
        if self.name:
            r = self.name + r
        return r

    __getattr__ = dict.__getitem__

    def m2io_repr(self):
        r = self.name + ' ' if self.name else ''
        r += '{\n'
        for name, value in self.items():
            if not isinstance(value, Table):
                r += name + '\n'
        r += '\n:::\n'
        for name, value in self.items():
            if not isinstance(value, Table):
                r += m2io_repr(value, name) + '\n'
        for name, value in self.items():
            if isinstance(value, Table):
                r += value.m2io_repr(name)
        r += '}\n'
        return r

def quotedsplit(handle):
    '''
    regex based string split. Handles double quotes and escapes inside
    double quotes.
    '''
    if not isinstance(handle, basestring):
        handle = handle.read()

    quoted_re = re.compile(r'"((?:\\.|[^"])*)"')
    comment_re = re.compile(r'\s#[^#]*#')
    escaped_re_sub = re.compile(r'\\(.)').sub

    reduced = quoted_re.sub('"', handle)
    reduced = comment_re.sub(' ', reduced)
    pos = 0

    for tok in reduced.split():
        if tok == '"':
            match = quoted_re.search(handle, pos)
            pos = match.end()
            tok = escaped_re_sub(r'\1', match.group(1))
        yield tok

class M2IOReader(object):
    '''
    Opens an M2IO file for reading and provides an iterator of all blocks.
    '''
    def __init__(self, filename):
        self.tokenizer = quotedsplit(open(filename))

    def __iter__(self):
        return self

    def aslist(self):
        '''
        Returns a list-like object of the remaining blocks in the file.
        '''
        return M2IO(self)

    def next(self):
        '''
        Reads and returns the next block.
        '''
        return self._next()

    def find_block(self, name):
        '''
        Find the next block with matching name.
        '''
        for block in self:
            if block.name == name:
                return block
        return None

    def _next(self, name=None):
        tok = self.tokenizer.next()
        if name is None and tok != '{':
            name = tok
            tok = self.tokenizer.next()
        assert tok == '{', tok + ' != {'
        items = zip(self._read_keys(), self.tokenizer)
        for tok in self.tokenizer:
            if tok == '}':
                break
            m = re_m.match(tok)
            if m is not None:
                items.append((m.group(1),
                    self._read_table(int(m.group(2)))))
            else:
                items.append((tok, self._next(tok)))
        return Block(name, items)

    def _read_keys(self):
        keys = []
        for tok in self.tokenizer:
            if tok == ':::':
                return keys
            keys.append(tok)
        raise ValueError

    def _read_table(self, nrows):
        tok = self.tokenizer.next()
        assert tok == '{', tok + ' != {'
        keys = self._read_keys()
        table = Table(keys)
        for tok in self.tokenizer:
            if tok == ':::':
                break
            table.add_row(self.tokenizer)
        tok = self.tokenizer.next()
        assert tok == '}', tok + ' != }'
        assert len(table) == nrows, 'len(table) != nrows'
        return table

if __name__ == '__main__':
    import sys
    import time

    for filename in sys.argv[1:]:
        start = time.time()

        try:
            M = M2IOReader(filename).aslist()
        except BaseException as e:
            print "failed", filename, e
        else:
            elapsed = time.time() - start
            print "succeded %5.2fs" % elapsed, filename, len(M)
