from util import Path

class Node (object):

	def __init__ (self, parent, path, exp, **kw):
		self.parent = parent
		self.suffix_link = None
		self.path = path
		self.name = kw.get('name', '')
		self.expGroup=exp[str_id]

	def __len__ (self):
		return len (self.path)

	def is_leaf (self): 
		return False

	def is_internal (self):
		return False

	def to_dot (self, a):
		if self.suffix_link is not None:
			a.append('"%s" -> "%s" [color=blue; constraint=false];\n' % (str (self), str (self.suffix_link)))


class Leaf (Node):
	def __init__ (self, parent, str_id, path, exp, **kw):
		super().__init__ (parent, path, **kw)
		self.str_id = str_id
		self.expGroup=exp[str_id]

	def __str__ (self):\
		return ("%s" % (self.name or str(self.path)) + ' %s:%d' % (self.str_id, self.path.start + 1))

	def is_leaf (self):
		return True

	def order (self, M):
		if self.is_leaf ():
			if self.path._end == Path.inf:
				self.path._end = M
		return

	def to_dot (self, a):
		a.append ('"%s" [color=green];\n' % str (self))
		super().to_dot (a)


class Internal (Node):
	def __init__ (self, parent, path,exp, **kw):
		super().__init__ (parent, path, **kw)
		self.children = {}
		self.expGroup=exp[str_id]

	def __str__ (self):\
		return "%s" % (self.name or str (self.path))

	def is_internal (self):
		return True

	def order (self, M):
		if self.is_leaf ():
			if self.path._end == Path.inf:
				self.path._end = M
		for node in self.children.values():
			node.order(M)
		return

	def split_edge (self, new_len, child):
		p1 = self.path
		p2 = child.path
		assert len (p1) < new_len < len (p2), "split length %d->%d->%d" % (
			len (p1), new_len, len (p2))
		edge_start = p2.start + len (p1)
		edge_end   = p2.start + new_len
		new_node = Internal (self, Path (p2.S, p2.start, edge_end))
		

		self.children[p2.S[edge_start]] = new_node     # substitute new node
		new_node.children [p2.S[edge_end  ]] = child
		child.parent = new_node

		return new_node

	def to_dot (self, a):
		a.append ('"%s" [color=red];\n' % str (self))
		super ().to_dot (a)
		for key in sorted (self.children):
			child = self.children[key]
			a.append ('"%s" -> "%s" [label="%s"];\n' % (str (self), str (child), str (key)))
			child.to_dot (a)
