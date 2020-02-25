import collections
import itertools
from nodes import *
from util import *

class Tree (object): 

	def __init__ (self, d, exp):
		""" 
		Given a dictionary of haplotype_id: list (each list corresponding to a haplotype),
		builds a generalized suffix tree using Ukkonen's algorithm.
		"""
		self.root = Internal (None, Path (tuple(), 0, 0), exp=self.exp, name = 'root')
		self.exp=exp
		self.leaves=[]

		for hap_id, haplotype in d.items ():
			self.id   = hap_id
			self.path = Path.construct(itertools.chain(haplotype, [UniqueEndChar(hap_id)]))
			self.build()


	def to_dot (self):
		"""
		"Draws" the tree in .dot format, to be used for debugging with test cases
		"""
		dot = []
		dot.append ('strict digraph G {\n')
		self.root.to_dot (dot)
		dot.append ('}\n')
		return ''.join (dot)

	def transition (self, node, k):

		if node is self.aux:
			return self.root, Path (self.path.S, 0, 1)

		node_prime = node.children[self.path[k]]
		path    = node_prime.path

		return node_prime, Path (path.S, path.start + len(node), path.end)


	def test_and_split (self, node, path, t):

		if path:
			node_prime, path_prime = self.transition(node, path.start)
			if t == path_prime[len (path)]:
				return True, node
			split_depth = len(node) + len (path)
			r = node.split_edge (split_depth, node_prime)
			return False, r
		else:
			if node is self.aux:
				return True, node
			if t in node.children:
				return True, node
			return False, node


	def canonize (self, node, path):
		if not path:
			return node, path

		node_prime, path_prime = self.transition (node, path.start)
		while len (path_prime) <= len (path):
			path.start += len (path_prime)
			node = node_prime

			if path:
				node_prime, path_prime = self.transition (node, path.start)

		return node, path


	def update (self, node, path):
		t_i = self.path[path.end-1]

		act_path = Path (path.S, path.start, path.end - 1)

		oldr = self.root
		is_end_point, r = self.test_and_split (node, act_path, t_i)
		while not is_end_point:
			start = path.end-1 - len (r)
			r_prime = Leaf (r, self.id, Path (self.path.S, start, Path.inf), exp=self.exp)
			if r_prime.is_leaf()==True:
				self.leaves.append(r_prime)
				print(r_prime.str_id)
			r.children[t_i] = r_prime

			if oldr is not self.root:
				oldr.suffix_link = r
			oldr = r

			node, act_path = self.canonize (node.suffix_link, act_path)
			is_end_point, r = self.test_and_split (node, act_path, t_i)

		if oldr is not self.root:
			oldr.suffix_link = node

		path.start = act_path.start
		return node, path


	def build (self):
		aux = Internal (None, Path (tuple(), 0, 0), exp=self.exp, name = 'aux')
		self.root.suffix_link = aux
		self.aux              = aux

		node    = self.root
		path = Path (self.path, 0, 1)
		len_ = len (self.path)

		while True:
			Path.e  = path.end
			node, path = self.update(node, path)
			node, path = self.canonize(node, path)

			if path.end == len_:
				break

			path.end += 1

		self.root.order(path.end)
		self.root.parent      = None
		self.root.suffix_link = None
		self.aux              = None

