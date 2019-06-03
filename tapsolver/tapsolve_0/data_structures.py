import networkx as nx
from networkx import Graph, DiGraph
import networkx.algorithms.isomorphism as iso
from networkx.readwrite import json_graph

from imolecule.notebook import generate
from imolecule.format_converter import convert
import imolecule
from ase import Atoms

def unique_name(namein, name_list):
    i = 1
    name0 = namein
    while namein in name_list:
        namein = name0 + '(' + str(i) + ')'
        i += 1
    return namein, i-1

class MolGraph(Graph):
    def __init__(self, data_dict=None):
        Graph.__init__(self)
        if data_dict is not None:
            self.from_dict(data_dict)
        self.duplicity = 0

    ### Read/Write functions ###
    def to_dict(self):
        """Construct molecule dict with atoms/bonds from molecule graph"""
        mol = {}
        atoms = []
        bonds = []
        nodes = []
        for n in self.nodes():
            if self.node[n]:
                atoms.append(self.node[n].copy())
                nodes.append(n)

        for e in self.edges():
            idxs = [nodes.index(e[0]), nodes.index(e[1])]
            bonds.append({'atoms': idxs, 'order': 1})
        mol['atoms'] = atoms
        mol['bonds'] = bonds

        return mol

    def to_atoms_object(self):
        """Construct molecule dict with atoms/bonds from molecule graph"""
        mol = {}
        atoms = []
        bonds = []
        nodes = []
        elements = []
        atom_coords = []
        for n in self.nodes():
            if self.node[n]:
                nodey = self.node[n].copy()
                atom_coords.append(nodey['location'])
                elements.append(nodey['element'])
                molecule = Atoms(elements,positions=atom_coords)
            return molecule

    def from_dict(self, moldict):#self, moldict):
        """Construct molecule graph from molecule dict with atoms/bonds"""
        bonds = moldict['bonds']
        atoms = moldict['atoms']
        nodes = range(len(moldict['atoms']))
        node_attr = []
        for ai in atoms:
            node = ai['element'] + str(ai['location'])
            attr = ai
            node_attr.append((node, attr))

        self.add_nodes_from(node_attr)

        for b in bonds:
            idxs = b['atoms']
            node_1 = node_attr[idxs[0]][0]
            node_2 = node_attr[idxs[1]][0]
            self.add_edge(node_1, node_2)

    def generate(self, data, data_format):
        """Wrapper around imolecule generate function"""
        moldict = eval(generate(data, data_format))
        self.from_dict(moldict)
        self.graph['chem_name'] = self.__altname__()
        #self.graph['chem_name'] = self.__altname__

        moldict.pop('molecular_weight')
        moldict.pop('formula')
        
        return self#,moldict

    def bonds(self):
        """Return bond names"""
        edges = self.edges()
        edge_names = []
        for a,b in edges:
            element_a = self.node[a]['element']
            element_b = self.node[b]['element']
            edge_name = '-'.join(sorted([element_a,element_b]))
            edge_names.append(edge_name)
        return edge_names

    def composition(self,composition=None):
        """Return composition dictionary.
        If composition dictionary is supplied
        then compute cumulatively."""
        if composition is None:
            composition = {}

        for n in self.nodes():
            e = self.node[n]['element']
            if e not in composition:
                composition[e] = 1
            else:
                composition[e] += 1
        return composition



    ### Data structure functions ###
    def __str__(self):
        """Shorthand linear text representation of the molecule"""
        txt = ''
        elements = sorted([(self.node[i]['element'],i) for i in self.nodes()])
        node_0 = elements[0][1]
        for i in nx.dfs_postorder_nodes(self, node_0):
            txt += self.node[i]['element']
        if self.duplicity:
            txt += '('+str(self.duplicity)+')'
        return txt

    #def check_oxygen(self,node,list,element):


    #ACY-in
    def __altname__(self):
        list_of_atoms = []
        iter_step = self.composition()
        if 'C' in iter_step:
            for n in self.nodes():
                if self.node[n]['element'] == 'C':
                    list_of_atoms.append('C')
                    H_nodes = 0
                    O_nodes = 0
                    OH_nodes= 0
                    N_nodes = 0
                    NO_nodes = 0
                    NH_nodes = 0
                    for a,b in self.edges(n):
                        if self.node[b]['element'] == 'H':
                            H_nodes += 1
                        if self.node[b]['element'] == 'O':
                            
                            #########################################
                            #####Need to correct this logic##########
                            #########################################
                            old_H_val = OH_nodes
                            for c,d in self.edges(b):
                                if self.node[d]['element'] == 'H':
                                    OH_nodes += 1
                            if OH_nodes-old_H_val == 0:
                                O_nodes += 1
                        if self.node[b]['element'] == 'N':
                            for c,d in self.edges(b):
                                if self.node[d]['element'] == 'H':
                                    NH_nodes += 1  
                                if self.node[d]['element'] == 'O':
                                    NO_nodes += 1
                                if NO_nodes and NH_nodes == 0:
                                    N_nodes += 1
                        #Need way to add check if a hydrogen atom comes after the oxygen
                        #list_of_atoms.append(self.node[b]['element'])
                    if H_nodes != 0:
                        if H_nodes > 1:
                            list_of_atoms.append('H'+str(H_nodes))
                        else:
                            list_of_atoms.append('H')
                    if O_nodes != 0:
                        if O_nodes > 1:
                            list_of_atoms.append('O'+str(O_nodes))
                        else:
                            list_of_atoms.append('O')
                    if OH_nodes != 0:
                        if OH_nodes > 1:
                            list_of_atoms.append('(OH)'+str(OH_nodes))
                        else:
                            list_of_atoms.append('OH')
                    if N_nodes != 0:
                        if N_nodes > 1:
                            list_of_atoms.append('N'+str(N_nodes))
                        else:
                            list_of_atoms.append('N')
                    if NO_nodes !=0:
                        if NO_nodes > 1:
                            list_of_atoms.append('NO'+str(NO_nodes))
                        else:
                            list_of_atoms.append('NO')
                    if NH_nodes !=0:
                        if NH_nodes > 1:
                            list_of_atoms.append('NH'+str(NH_nodes))
                        else:
                            list_of_atoms.append('NH')

        elif 'N' in iter_step:
            for n in self.nodes():
                if self.node[n]['element'] == 'N':
                    list_of_atoms.append('N')
                    H_nodes = 0
                    O_nodes = 0
                    OH_nodes= 0
                    N_nodes = 0
                    NO_nodes = 0
                    NH_nodes = 0
                    for a,b in self.edges(n):
                        if self.node[b]['element'] == 'H':
                            H_nodes += 1
                        if self.node[b]['element'] == 'O':
                            
                            #########################################
                            #####Need to correct this logic##########
                            #########################################
                            old_H_val = OH_nodes
                            for c,d in self.edges(b):
                                if self.node[d]['element'] == 'H':
                                    OH_nodes += 1
                            if OH_nodes-old_H_val == 0:
                                O_nodes += 1
                    if H_nodes != 0:
                        if H_nodes > 1:
                            list_of_atoms.append('H'+str(H_nodes))
                        else:
                            list_of_atoms.append('H')
                    if O_nodes != 0:
                        if O_nodes > 1:
                            list_of_atoms.append('O'+str(O_nodes))
                        else:
                            list_of_atoms.append('O')
                    if OH_nodes != 0:
                        if OH_nodes > 1:
                            list_of_atoms.append('(OH)'+str(OH_nodes))
                        else:
                            list_of_atoms.append('OH')
        elif 'O' in iter_step:
            step_i = 0
            for n in self.nodes():
                if self.node[n]['element'] == 'O':
                    if step_i == 1:
                        break
                    step_i += 1
                    OH_nodes= 0
                    H_nodes = 0
                    O_nodes = 0
                    for a,b in self.edges(n):
                        if self.node[b]['element'] == 'H':
                            H_nodes += 1
                        if self.node[b]['element'] == 'O':
                            
                            #########################################
                            #####Need to correct this logic##########
                            #########################################
                            old_H_val = OH_nodes
                            for c,d in self.edges(b):
                                if self.node[d]['element'] == 'H':
                                    OH_nodes += 1
                            if OH_nodes-old_H_val == 0:
                                O_nodes += 1
                    if H_nodes != 0:
                        if H_nodes > 1:
                            list_of_atoms.append('H2O')
                        else:
                            list_of_atoms.append('OH')
                    if O_nodes != 0:
                        list_of_atoms.append('O2')
                    if OH_nodes != 0:
                        if OH_nodes > 1:
                            list_of_atoms.append('OOH')
                    if len(list_of_atoms) == 0:
                        list_of_atoms.append('O')
        else:
            if self.number_of_nodes() > 1:
                list_of_atoms.append('H2')
            else:
               list_of_atoms.append('H')    
        test_print = ''.join(list_of_atoms)
        return test_print
    #ACY-out

    def __repr__(self):
        """Comprehensive representation of object. eval(repr(x)) == x."""
        return str('MolGraph({})'.format(self.to_dict()))

    def __eq__(self, other):
        """Check equivalencey between two molecules. Note that this is
        topological equivalency (e.g. atom/bonds) as defined by graph
        isomorphism. Atoms are considered equivalent if the element and
        charge are the same."""
        nm = iso.categorical_node_match('element', 'charge')
        return nx.is_isomorphic(self, other, node_match=nm)

class RxnGraph(DiGraph):
    def __init__(self,*args,**kwargs):
        DiGraph.__init__(self,*args,**kwargs)
        mol_defaults = {'size': 1,
                        'color':'black',
                        'label':'{name}',
                        'description':'{name}',
                        }
        rxn_defaults = {'width': 0.1,
                        'height':0.1,
                        'fixed_size':'true',
                        'label':'',
                        'description':'{name}',
                        'penwidth':3,
                        }
        self.graphviz_format = {}
        self.graphviz_format['molecule_defaults'] = mol_defaults
        self.graphviz_format['reaction_defaults'] = rxn_defaults

    def add_molecule(self, node, attrs={}):
        """Add molecule to reaction graph.
        Ensure that it is not a duplicate by checking
        against existing molecule graphs"""

        if 'graph' not in attrs:
            raise AttributeError('Molecules must have associated graph.')

        if not hasattr(attrs['graph'],'is_directed'):
            raise AttributeError('Molecule graph must be graph object.')

        molecules = [n for n in self.nodes() if self.node[n]['type'] == 'molecule']

        graph_matches = [m for m in molecules if self.node[m]['graph'] == attrs['graph']]

        name_matches = [m for m in molecules if m == node]

        if not graph_matches and not name_matches: 
            #molecule name and graph are unique
            #DiGraph.add_node(self,node,attrs)
            DiGraph.add_node(self, node, attrs)
            return node

        elif not graph_matches and len(name_matches) == 1:
            #molecule is unique, but name is taken
            new_name, index = unique_name(node,self.nodes())
            attrs['graph'].duplicity = index
            DiGraph.add_node(self,new_name,attrs)
            return new_name

        elif not name_matches and len(graph_matches) == 1:
            #molecule is in graph, but has a different name
            node_name = graph_matches[0]
            return node_name

        elif len(name_matches) == 1 and len(graph_matches) == 1:
            if self.node[node]['graph'] == attrs['graph']:
                #node is already in the graph
                return node
            else:
                #molecule is in the graph but there is a name conflict.
                #return the proper name.
                node_name = graph_matches[0]
                return node_name

        elif len(graph_matches) > 1:
            raise ValueError('Molecule {} already has a duplicate in the graph. Always use add_molecule rather than add_node to avoid this error.'.format(node))

        else:
            raise ValueError('Unknown error')

    ### These read/write functions need cleanup ###
    def from_rxn_list(self, rxn_list):
        react_num = 0
        for rxn in rxn_list:
            reactants, products = rxn
            reaction_number = []
            r_names = []
            p_names = []
            for r in reactants + products:
                assert hasattr(r, 'is_directed')  # ensure that r is MolGraph
                attrs = {'type': 'molecule', 'graph': r}
                if r in reactants:
                    attrs['molecule_type'] = 'reactant'
                    attrs['chem_name'] = r.__altname__()
                    r_name = self.add_molecule(str(r.graph['chem_name']), attrs)
                    r_names.append(r_name)
                if r in products:
                    attrs['molecule_type'] = 'product'
                    r.graph['chem_name'] = r.__altname__()
                    r_name = self.add_molecule(str(r.graph['chem_name']), attrs)
                    p_names.append(r_name)
            #Add here
            rxn_name = '+'.join([str(r)+'*' for r in r_names]) #!!!!!!!!!!!!!!!!
            rxn_name += ' <-> '
            rxn_name += ' + '.join([str(p)+'*' for p in p_names])

            if rxn_name in self.nodes():
                raise ValueError('Duplicate reaction detected. Using add_molecule instead of add_node to add molecules should avoid this error.')

            self.add_node(rxn_name, attr_dict={'type': 'reaction'})

            for r_name in r_names:
                self.add_edge(r_name, rxn_name)
            for r_name in p_names:
                self.add_edge(rxn_name, r_name)

            #print(rxn_name)
        

    def to_rxn_list(self):
        rxns = [r for r in self.nodes() if self.node[r]['type'] == 'reaction']
        all_rxns = []
        for rxn in rxns:
            prods = [self.node[p]['chem_name'] for p in self.successors(rxn)]
            reacts = [self.node[p]['chem_name'] for p in self.predecessors(rxn)]
            all_rxns.append([reacts, prods])
        return all_rxns

    def format_graphviz(self):
        """Format for visualization with graphviz"""
        for n in self.nodes():
            attrs = {}
            if self.node[n]['type'] == 'reaction':
                name = n
                defaults = self.graphviz_format['reaction_defaults']
                attrs = defaults.copy()
                for key in attrs:
                    if hasattr(attrs[key],'format') and '{' in attrs[key]:
                        attrs[key] = attrs[key].format(**locals())
                node_type = self.node[n]['reaction_type'] 

            elif self.node[n]['type'] == 'molecule':
                name = n
                defaults = self.graphviz_format['molecule_defaults']
                attrs = defaults.copy()
                for key in attrs:
                    if hasattr(attrs[key],'format') and '{' in attrs[key]:
                        attrs[key] = attrs[key].format(**locals())
                node_type = self.node[n]['molecule_type'] 

            if node_type in self.graphviz_format:
                attrs.update(self.graphviz_format[node_type])

            self.node[n].update(attrs)

            if self.node[n]['type'] == 'reaction':
                edges = self.in_edges([n]) + self.out_edges([n])
                for u,v in edges:
                    self[u][v] = attrs

    def to_agraph(self):
        """Return pygraphviz agraph string for 2D drawing"""

        return str(nx.drawing.nx_agraph.to_agraph(self))

    def to_jgraph(self):
        """Return json object for jgraph rendering"""
        jg = {'nodes':{}, 'edges':[]}
        node_name_dict = {}
        for n in self.nodes():
            attrs = {}
            if self.node[n]['type'] == 'reaction':
                attrs['size'] = 0.25
                attrs['color'] = '0x000000'
                name = n
                attrs['label'] = name
                attrs['.label'] = name
                attrs['description'] = name
                node_name_dict[n] = name
            elif self.node[n]['type'] == 'molecule':
                attrs['size'] = 1
                attrs['color'] = '0x0000ff'
                name = n
                i = 1
                while name in jg['nodes']:
                    name = n+'('+str(i)+')'

                attrs['label'] = name
                attrs['.label'] = name
                attrs['description'] = name
                node_name_dict[n] = name

            jg['nodes'][name] = attrs

        for e in self.edges():
            source, target = e
            sn = node_name_dict[source]
            tn = node_name_dict[target]
            jg['edges'].append({'source':sn, 'target':tn})

        return jg

    def reverse_rxn(self, n, remove=False):
        if not self.node[n]['type'] == 'reaction':
            print('Only reaction nodes can be reversed.')
            return None

        reacts, prods = self.get_reactants_products(n)
        react_edges = [(str(r), n) for r in reacts]
        prod_edges = [(n,str(p)) for p in prods]
        relevant_edges = react_edges+prod_edges
        rxn_class = self.node[n].get('reaction_type',None)

        if remove:
            self.remove_edges_from(relevant_edges)
            rxn_name = '+'.join([str(r) for r in prods])
            rxn_name += '->'
            rxn_name += '+'.join([str(p) for p in reacts])
            if rxn_class:
                if 'coupling' in rxn_class:
                    rxn_class = rxn_class.replace('coupling','scission')
                if 'scission' in rxn_class:
                    rxn_class = rxn_class.replace('scission','coupling')
                self.node[n]['reaction_type'] = rxn_class
        else:
            rxn_name = '+'.join([str(r) for r in prods])
            rxn_name += '<->'
            rxn_name += '+'.join([str(p) for p in reacts])
            if rxn_class:
                if 'coupling' in rxn_class:
                    rxn_class += '+scission'
                if 'scission' in rxn_class:
                    rxn_class += '+coupling'
                self.node[n]['reaction_type'] = rxn_class

        mapping = {n:rxn_name}
        nx.relabel.relabel_nodes(self,mapping,copy=False)
        
        reversed_edges = [(str(r),rxn_name) for r in prods]
        reversed_edges += [(rxn_name,str(p)) for p in reacts]
        self.add_edges_from(reversed_edges)

    def get_reactants_products(self,n):
        if not self.node[n]['type'] == 'reaction':
            print('Only reaction nodes have reactants/products.')
            return None, None

        prods = [self.node[p]['graph'] for p in self.successors(n)]
        reacts = [self.node[p]['graph'] for p in self.predecessors(n)]
        return reacts, prods

    def classify_rxns(self,verbose=False):
        rxns = [r for r in self.nodes() if (self.node[r]['type'] == 'reaction' and self.node[r].get('reaction_type',None) == None)] #only unclassified reactions
        for rxn in rxns:
            reacts, prods = self.get_reactants_products(rxn)
            if verbose:
                print(rxn)
                print('Reactants: '+str([str(r) for r in reacts]))
                print('Products: '+str([str(r) for r in prods]))

            if len(reacts) == len(prods) and '<->' in rxn: # reversible reaction
                print('Warning: Classification for reversible reactions relies on reaction node names and may not be robust. Reactions should be classified before they are made reversible.')
                if '(' in rxn:
                    rxn_name = rxn.split('(')[0]
                else:
                    rxn_name = rxn
                react_names,prod_names = rxn_name.split('<->')
                react_names = react_names.split('+')
                prod_names = prod_names.split('+')
                reacts = [r for r in reacts if str(r) in react_names]
                prods = [p for p in prods if str(p) in prod_names]

            react_edges = []
            prod_edges = []

            prod_composition = {} #use compositions to check for symmetric scissions
            react_composition = {}
            for r in reacts:
                print(r)
                print(type(r))
                graph = self.node[r]['graph'] #self.node[str(r)]['graph']
                print(graph)
                react_composition = graph.composition(react_composition)
                react_edges += graph.bonds()

            for p in prods:
                graph = self.node[str(p)]['graph']
                prod_composition = graph.composition(prod_composition)
                prod_edges += graph.bonds()

            symmetric_scission = True
            for element in react_composition:
                if (react_composition[element] != prod_composition.get(element,0)*2):
                    symmetric_scission = False
            if symmetric_scission:
                prod_edges *= 2

            symmetric_coupling = True
            for element in prod_composition:
                if (prod_composition[element] != react_composition.get(element,0)*2):
                    symmetric_coupling = False
            if symmetric_coupling:
                react_edges *= 2

            if len(prod_edges) > len(react_edges):
                diff = [e for e in prod_edges]
                for e in react_edges:
                    if e in diff:
                        diff.remove(e)
                rxn_class = 'coupling'

            elif len(prod_edges) < len(react_edges):
                diff = [e for e in react_edges]
                for e in prod_edges:
                    if e in diff:
                        diff.remove(e)
                rxn_class = 'scission'
            else:
                rxn_class = 'isomerization'
                diff = []

            if len(diff) > 1:
                print('Warning: Could not reliably classify {}. Assuming {}.'.format(rxn,rxn_class))
                print(react_edges)
                print(prod_edges)
                print(diff)
            elif len(diff) == 1:
                rxn_class = diff[0]+'_'+rxn_class

            self.node[rxn]['reaction_type'] = rxn_class

    @staticmethod
    def node_matcher(n1, n2):
        """Define whether or not nodes are equivalent. Static because
        it does not depend on the properties of the entire graph"""
        if n1['type'] != n2['type']:
            # reactions are not the same as molecules
            return False
        elif n1['type'] == n2['type'] == 'molecule':
            # molecules are only equivalent if their graphs are.
            # see MolGraph.__eq__
            return n1['graph'] == n2['graph']
        elif n1['type'] == n2['type'] == 'reaction':
            # treat all reaction nodes as equivalent
            return True

    def __str__(self):
        """Shorthand text representation of the reactions"""
        rxns = []
        for n in self.nodes():
            if self.node[n]['type'] == 'reaction':
                rxns.append(str(n))
        txt = '\n'.join(rxns)
        return txt

    def __repr__(self):
        data = json_graph.node_link_data(self)
        return repr(data)

    def __eq__(self, other):
        """Check equivalencey between two reaction networks.
        Molecules are considered equivalent based on topology."""
        return nx.is_isomorphic(self, other, node_match=self.node_matcher)