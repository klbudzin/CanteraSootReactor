"""
 Cantera .cti input file processor

 The functions and classes in this module process Cantera .cti input
 files and produce CTML files. It can be imported as a module, or used
 as a script.

 script usage:

 python ctml_writer.py infile.cti

 This will produce CTML file 'infile.xml'

"""
import string

class CTI_Error:
    """Exception raised if an error is encountered while
    parsing the input file."""
    def __init__(self, msg):
        print '\n\n***** Error parsing input file *****\n\n'
        print msg
        print




indent = ['',
          ' ',
          '  ',
          '   ',
          '    ',
          '     ',
          '      ',
          '       ',
          '        ',
          '          ',
          '           ',
          '            ',
          '             ',
          '              ',
          '               ',
          '                ']

#-----------------------------------------------------

class XMLnode:
    
    """This is a minimal class to allow easy creation of an XML tree
    from Python. It can write XML, but cannot read it."""
    
    def __init__(self, name="--", value = ""):
        
        """Create a new node. Usually this only needs to be explicitly
        called to create the root element. Method addChild calls this
        constructor to create the new child node."""

        self._name = name

        # convert 'value' to a string if it is not already, and
        # strip leading whitespace        
        if type(value) <> types.StringType:
            self._value = string.lstrip(`value`)
        else:
            self._value = string.lstrip(value)
            
        self._attribs = {}    # dictionary of attributes
        self._children = []   # list of child nodes
        self._childmap = {}   # dictionary of child nodes

        
    def name(self):
        """The tag name of the node."""
        return self._name

    def nChildren(self):
        """Number of child elements."""
        return len(self._children)

    def addChild(self, name, value=""):
        """Add a child with tag 'name', and set its value if the value
        parameter is supplied."""
        
        # create a new node for the child
        c = XMLnode(name = name, value = value)
        
        # add it to the list of children, and to the dictionary
        # of children
        self._children.append(c)
        self._childmap[name] = c
        return c

    def addComment(self, comment):
        """Add a comment."""
        self.addChild(name = '_comment_', value = comment)
        
    def value(self):
        """A string containing the element value."""
        return self._value

    def child(self, name=""):
        """The child node with specified name."""
        return self._childmap[name]

    def __getitem__(self, key):
        """Get an attribute using the syntax node[key]"""
        return self._attribs[key]

    def __setitem__(self, key, value):
        """Set a new attribute using the syntax node[key] = value."""
        self._attribs[key] = value    

    def __call__(self):
        """Allows getting the value using the syntax 'node()'"""
        return self._value

    def write(self, file):
        """Write out the XML tree to a file."""
        f = open(file,'w')
        f.write('<?xml version="1.0"?>\n')
        self._write(f, 0)
        f.write('\n')
        
    def _write(self, f, level = 0):

        """Internal method used to write the XML representation of
        each node."""
        if self._name == "": return
        
        indnt = indent[level]

        # handle comments
        if self._name == '_comment_':
            f.write('\n'+indnt+'<!--')
            if len(self._value) > 0:
                if self._value[0] <> ' ':
                    self._value = ' '+self._value
                if self._value[-1] <> ' ':
                    self._value = self._value+' '
            f.write(self._value+'-->')
            return

        # write the opening tag and attributes
        f.write(indnt + '<' + self._name)
        for a in self._attribs.keys():
            f.write(' '+a+'="'+self._attribs[a]+'"')
        if (self._value == "" and self.nChildren() == 0):
            f.write('/>')
        else:
            f.write('>')
            if self._value <> "":
                vv = string.lstrip(self._value)
                ieol = vv.find('\n')
                if ieol >= 0:
                    while 1 > 0:
                        ieol = vv.find('\n')
                        if ieol >= 0:
                            f.write('\n  '+indnt+vv[:ieol])
                            vv = string.lstrip(vv[ieol+1:])
                        else:
                            f.write('\n  '+indnt+vv)
                            break
                else:
                    f.write(self._value)

            for c in self._children:
                f.write('\n')
                c._write(f, level + 2)
            if (self.nChildren() > 0):
                f.write('\n'+indnt)
            f.write('</'+self._name+'>')

#--------------------------------------------------

# constants that can be used in .cti files
OneAtm = 1.01325e5
OneBar = 1.0e5


import types, math, copy

# default units
_ulen = 'm' 
_umol = 'kmol'
_umass = 'kg'
_utime = 's'
_ue = 'J/kmol'
_uenergy = 'J'
_upres = 'Pa'

# used to convert reaction pre-exponentials
_length = {'cm':0.01, 'm':1.0, 'mm':0.001}
_moles = {'kmol':1.0, 'mol':0.001, 'molec':1.0/6.023e26}
_time = {'s':1.0, 'min':60.0, 'hr':3600.0}

# default std state pressure
_pref = 1.0e5    # 1 bar

_name = 'noname'

# these lists store top-level entries
_species = []
_speciesnames = []
_phases = []
_reactions = []
_atw = {}


_valsp = ''
_valrxn = ''


def validate(species = 'yes', reactions = 'yes'):
    global _valsp
    global _valrxn
    _valsp = species
    _valrxn = reactions
        
def isnum(a):
    """True if a is an integer or floating-point number."""
    if type(a) == types.IntType or type(a) == types.FloatType:
        return 1
    else:
        return 0

def is_local_species(name):
    """true if the species named 'name' is defined in this file"""
    if name in _speciesnames:
        return 1
    return 0

def dataset(nm):
    "Set the dataset name. Invoke this to change the name of the xml file."
    global _name
    _name = nm

def standard_pressure(p0):
    """Set the default standard-state pressure."""
    global _pref
    _pref = p0

def units(length = '', quantity = '', mass = '', time = '',
          act_energy = '', energy = '', pressure = ''):
    """set the default units."""
    global _ulen, _umol, _ue, _utime, _umass, _uenergy, _upres
    if length: _ulen = length
    if quantity: _umol = quantity
    if act_energy: _ue = act_energy
    if time: _utime = time
    if mass: _umass = mass
    if energy: _uenergy = energy
    if pressure: _upres = pressure

def ufmt(base, n):
    """return a string representing a unit to a power n."""
    if n == 0: return ''
    if n == 1: return '-'+base
    if n == -1: return '/'+base
    if n > 0: return '-'+base+`n`
    if n < 0: return '/'+base+`-n`    
    
def write():
    """write the CTML file."""
    x = XMLnode("ctml")
    v = x.addChild("validate")
    v["species"] = _valsp
    v["reactions"] = _valrxn
    
    for ph in _phases:
        ph.build(x)
    s = species_set(name = _name, species = _species)
    s.build(x)
    
    r = x.addChild('reactionData')
    r['id'] = 'reaction_data'
    for rx in _reactions:
        rx.build(r)
    
    if _name <> 'noname':
        x.write(_name+'.xml')
    else:
        print x

            
def addFloat(x, nm, val, fmt='', defunits=''):
    """
    Add a child element to XML element x representing a
    floating-point number.
    """
    u = ''
    s = ''
    if isnum(val):
        fval = float(val)
        if fmt:
            s = fmt % fval
        else:
            s = `fval`
        xc = x.addChild(nm, s)
        if defunits:
            xc['units'] = defunits
    else:
        v = val[0]
        u = val[1]
        if fmt:
            s = fmt % v
        else:
            s = `v`
        xc = x.addChild(nm, s)
        xc['units'] = u

    
def getAtomicComp(atoms):
    if type(atoms) == types.DictType: return atoms
    a = atoms.replace(',',' ')
    toks = a.split()
    d = {}
    for t in toks:
        b = t.split(':')
        d[b[0]] = int(b[1])
    return d

def getReactionSpecies(s):
    toks = s.replace(' + ',' ').split()
    d = {}
    n = 1
    for t in toks:
        if t > '0' and t < '9':
            n = int(t)
        else:
            if d.has_key(t):
                d[t] += n
            else:
                d[t] = n
            n = 1
    return d

class species_set:
    def __init__(self, name = '', species = []):
        self._s = species
        self._name = name
        #self.type = SPECIES_SET
        
    def build(self, p):
        p.addComment('     species definitions     ')
        sd = p.addChild("speciesData")
        sd["id"] = "species_data"
        for s in self._s:
            #if s.type == SPECIES:
            s.build(sd)
            #else:
            #    raise 'wrong object type in species_set: '+s.__class__

            
class species:
    """A species."""
    
    def __init__(self,
                 name = 'missing name!',
                 atoms = '',
                 note = '',
                 thermo = None,
                 transport = None,
                 charge = -999):
        self._name = name
        self._atoms = getAtomicComp(atoms)

        self._comment = note
        
        if thermo:
            self._thermo = thermo
        else:
            self._thermo = const_cp()
            
        self._transport = transport
        chrg = 0
        self._charge = charge
        if self._atoms.has_key('E'):
            chrg = -self._atoms['E']
            if self._charge <> -999:
                if self._charge <> chrg:
                    raise CTI_Error('specified charge inconsistent with number of electrons')
            else:
                self._charge = chrg
                #        self.type = SPECIES
        
        global _species
        _species.append(self)
        global _speciesnames
        if name in _speciesnames:
            raise CTI_Error('species '+name+' multiply defined.')
        _speciesnames.append(name)
        

    def build(self, p):
        hdr = '    species '+self._name+'    '
        p.addComment(hdr)        
        s = p.addChild("species")
        s["name"] = self._name
        a = ''
        for e in self._atoms.keys():
            a += e+':'+`self._atoms[e]`+' '
        s.addChild("atomArray",a)
        if self._comment:
            s.addChild("note",self._comment)
        if self._charge <> -999:
            s.addChild("charge",self._charge)
        if self._thermo:
            t = s.addChild("thermo")
            if type(self._thermo) == types.InstanceType:
                self._thermo.build(t)
            else:
                nt = len(self._thermo)
                for n in range(nt):
                    self._thermo[n].build(t)
        if self._transport:
            t = s.addChild("transport")
            if type(self._transport) == types.InstanceType:
                self._transport.build(t)
            else:
                nt = len(self._transport)
                for n in range(nt):
                    self._transport[n].build(t)                    

class thermo:
    """Base class for species standard-state thermodynamic properties."""
    def _build(self, p):
        return p.addChild("thermo")

        
class NASA(thermo):
    """NASA polynomial parameterization."""
    
    def __init__(self, range = (0.0, 0.0), 
                 coeffs = [], p0 = -1.0):
        self._t = range
        self._pref = p0
        if len(coeffs) <> 7:
            raise CTI_Error('NASA coefficient list must have length = 7')
        self._coeffs = coeffs

        
    def build(self, t):
        n = t.addChild("NASA")
        n['Tmin'] = `self._t[0]`
        #n['Tmid'] = `self._t[1]`        
        n['Tmax'] = `self._t[1]`
        if self._pref <= 0.0:
            n['P0'] = `_pref`
        else:
            n['P0'] = `self._pref`
        str = ''
        for i in  range(4):
            str += '%17.9E, ' % self._coeffs[i]
        str += '\n'
        str += '%17.9E, %17.9E, %17.9E' % (self._coeffs[4],
                                           self._coeffs[5], self._coeffs[6])
        #if i > 0 and 3*((i+1)/3) == i: str += '\n'
        #str = str[:-2]
        u = n.addChild("floatArray", str)
        u["size"] = "7"
        u["name"] = "coeffs"


class Shomate(thermo):
    """Shomate polynomial parameterization."""
    
    def __init__(self, range = (0.0, 0.0), 
                 coeffs = [], p0 = -1.0):
        self._t = range
        self._pref = p0
        if len(coeffs) <> 7:
            raise CTI_Error('Shomate coefficient list must have length = 7')
        self._coeffs = coeffs

        
    def build(self, t):
        n = t.addChild("Shomate")
        n['Tmin'] = `self._t[0]`
        n['Tmax'] = `self._t[1]`
        if self._pref <= 0.0:
            n['P0'] = `_pref`
        else:
            n['P0'] = `self._pref`
        str = ''
        for i in  range(4):
            str += '%17.9E, ' % self._coeffs[i]
        str += '\n'
        str += '%17.9E, %17.9E, %17.9E' % (self._coeffs[4],
                                           self._coeffs[5], self._coeffs[6])
        u = n.addChild("floatArray", str)
        u["size"] = "7"
        u["name"] = "coeffs"
    
                 
class const_cp(thermo):
    """Constant specific heat."""
    
    def __init__(self, 
                 t0 = 298.15, cp0 = 0.0, h0 = 0.0, s0 = 0.0,
                 tmax = 5000.0, tmin = 100.0):
        self._t = [tmin, tmax]
        self._c = [t0, h0, s0, cp0]
        
    def build(self, t):
        #t = self._build(p)
        c = t.addChild('const_cp')
        if self._t[0] >= 0.0: c['Tmin'] = `self._t[0]`
        if self._t[1] >= 0.0: c['Tmax'] = `self._t[1]`
        energy_units = _uenergy+'/'+_umol
        addFloat(c,'t0',self._c[0], defunits = 'K')
        addFloat(c,'h0',self._c[1], defunits = energy_units)
        addFloat(c,'s0',self._c[2], defunits = energy_units+'/K')        
        addFloat(c,'cp0',self._c[3], defunits = energy_units+'/K')


class gas_transport:
    """Transport coefficients for ideal gas transport model."""

    def __init__(self, geom = 'nonlin',
                 diam = 0.0, well_depth = 0.0, dipole = 0.0,
                 polar = 0.0, rot_relax = 0.0):
        self._geom = geom
        self._diam = diam
        self._well_depth = well_depth
        self._dipole = dipole
        self._polar = polar
        self._rot_relax = rot_relax

    def build(self, t):
        #t = s.addChild("transport")
        t['model'] = 'gas_transport'
        #        t.addChild("geometry", self._geom)
        tg = t.addChild('string',self._geom)
        tg['title'] = 'geometry'
        addFloat(t, "LJ_welldepth", (self._well_depth, 'K'), '%8.3f')
        addFloat(t, "LJ_diameter", (self._diam, 'A'),'%8.3f')
        addFloat(t, "dipoleMoment", (self._dipole, 'Debye'),'%8.3f')        
        addFloat(t, "polarizability", (self._polar, 'A3'),'%8.3f')
        addFloat(t, "rotRelax", self._rot_relax,'%8.3f')        

        
class Arrhenius:
    def __init__(self,
                 A = 0.0,
                 n = 0.0,
                 E = 0.0,
                 coverage = [],
                 rate_type = ''):
        self._c = [A, n, E]
        self._type = rate_type
        
        if coverage:
            if type(coverage[0]) == types.StringType:
                self._cov = [coverage]
            else:
                self._cov = coverage
        else:
            self._cov = None

        
    def build(self, p, units_factor = 1.0,
              gas_species = [], name = '', rxn_phase = None):
        
        a = p.addChild('Arrhenius')
        if name: a['name'] = name

        # check for sticking probability
        if self._type:
            a['type'] = self._type
            if self._type == 'stick':
                ngas = len(gas_species)
                if ngas <> 1:
                    raise CTI_Error("""
Sticking probabilities can only be used for reactions with one gas-phase
reactant, but this reaction has """+`ngas`+': '+`gas_species`)
                else:
                    a['species'] = gas_species[0]
                    units_factor = 1.0
                
                    
        # if a pure number is entered for A, multiply by the conversion
        # factor to SI and write it to CTML as a pure number. Otherwise,
        # pass it as-is through to CTML with the unit string.
        if isnum(self._c[0]):
            addFloat(a,'A',self._c[0]*units_factor, fmt = '%14.6E')
        elif len(self._c[0]) == 2 and self._c[0][1] == '/site':
            addFloat(a,'A',self._c[0][0]/rxn_phase._sitedens,
                     fmt = '%14.6E')            
        else:
            addFloat(a,'A',self._c[0], fmt = '%14.6E')

        # The b coefficient should be dimensionless, so there is no
        # need to use 'addFloat'
        a.addChild('b',`self._c[1]`)

        # If a pure number is entered for the activation energy,
        # add the default units, otherwise use the supplied units.
        addFloat(a,'E', self._c[2], fmt = '%f', defunits = _ue)

        # for surface reactions, a coverage dependence may be specified.
        if self._cov:
            for cov in self._cov:
                c = a.addChild('coverage')
                c['species'] = cov[0]
                addFloat(c, 'a', cov[1], fmt = '%f')
                c.addChild('m', `cov[2]`)
                addFloat(c, 'e', cov[3], fmt = '%f', defunits = _ue)

def stick(A = 0.0, n = 0.0, E = 0.0, coverage = []):
    return Arrhenius(A = A, n = n, E = E, coverage = coverage, rate_type = 'stick')


def getPairs(s):
    toks = s.split()
    m = {}
    for t in toks:
        key, val = t.split(':')
        m[key] = float(val)
    return m
    
class reaction:
    def __init__(self,
                 equation = '',
                 kf = None,
                 id = '',
                 order = '',
                 options = []
                 ):

        self._id = id
        self._e = equation
        self._order = order
        
        if type(options) == types.StringType:
            self._options = [options]
        else:
            self._options = options
        global _reactions
        self._num = len(_reactions)+1
        r = ''
        p = ''
        for e in ['<=>', '=>', '=']:
            if self._e.find(e) >= 0:
                r, p = self._e.split(e)
                if e in ['<=>','=']: self.rev = 1
                else: self.rev = 0
                break
        self._r = getReactionSpecies(r)
        self._p = getReactionSpecies(p)

        self._rxnorder = copy.copy(self._r)
        if self._order:
            ord = getPairs(self._order)
            for o in ord.keys():
                if self._rxnorder.has_key(o):
                    self._rxnorder[o] = ord[o]
                else:
                    raise CTI_Error("order specified for non-reactant: "+o)
            
        self._kf = kf
        self._igspecies = []
        self._dims = [0]*4
        self._rxnphase = None
        self._type = ''
        _reactions.append(self)

        
    def build(self, p):
        if self._id:
            id = self._id
        else:
            if self._num < 10:
                nstr = '000'+`self._num`
            elif self._num < 100:
                nstr = '00'+`self._num`
            elif self._num < 1000:
                nstr = '0'+`self._num`                
            else:
                nstr = `self._num`
            id = nstr

        
        mdim = 0
        ldim = 0
        str = ''

        rxnph = []
        for s in self._r.keys():
            ns = self._rxnorder[s]
            nm = -999
            nl = -999

            str += s+':'+`self._r[s]`+' '
            mindim = 4
            for ph in _phases:
                if ph.has_species(s):
                    nm, nl = ph.conc_dim()
                    if ph.is_ideal_gas():
                        self._igspecies.append(s)
                    if not ph in rxnph:
                        rxnph.append(ph)
                        self._dims[ph._dim] += 1
                        if ph._dim < mindim:
                            self._rxnphase = ph
                            mindim = ph._dim
                    break
            if nm == -999:
                raise CTI_Error("species "+s+" not found")

            mdim += nm*ns
            ldim += nl*ns

        p.addComment("   reaction "+id+"    ")                
        r = p.addChild('reaction')
        r['id'] = id
        if self.rev:
            r['reversible'] = 'yes'
        else:
            r['reversible'] = 'no'
                
        noptions = len(self._options)
        for nss in range(noptions):
            s = self._options[nss]
            if s == 'duplicate':
                r['duplicate'] = 'yes'
            elif s == 'negative_A':
                r['negative_A'] = 'yes'

        ee = self._e.replace('<','[')
        ee = ee.replace('>',']')
        r.addChild('equation',ee)

        if self._order:
            for osp in self._rxnorder.keys():
                o = r.addChild('order',self._rxnorder[osp])
                o['species'] = osp
                

        # adjust the moles and length powers based on the dimensions of
        # the rate of progress (moles/length^2 or moles/length^3)
        if self._type == 'surface':
            mdim += -1
            ldim += 2
            p = self._dims[:3]
            if p[0] <> 0 or p[1] <> 0 or p[2] > 1:            
                raise CTI_Error(self._e +'\nA surface reaction may contain at most '+
                                   'one surface phase.')
        elif self._type == 'edge':
            mdim += -1
            ldim += 1
            p = self._dims[:2]
            if p[0] <> 0 or p[1] > 1:
                raise CTI_Error(self._e+'\nAn edge reaction may contain at most '+
                                   'one edge phase.')                
        else:
            mdim += -1
            ldim += 3

        # add the reaction type as an attribute if it has been specified.
        if self._type:
            r['type'] = self._type

        # The default rate coefficient type is Arrhenius. If the rate
        # coefficient has been specified as a sequence of three
        # numbers, then create a new Arrhenius instance for it;
        # otherwise, just use the supplied instance.
        nm = ''
        kfnode = r.addChild('rateCoeff')
        if self._type == '':
            self._kf = [self._kf]
        elif self._type == 'surface':
            self._kf = [self._kf]
        elif self._type == 'edge':
            self._kf = [self._kf]                        
        elif self._type == 'threeBody':
            self._kf = [self._kf]
            mdim += 1
            ldim -= 3

        if self._type == 'edge':
            if self._beta > 0:
                electro = kfnode.addChild('electrochem')
                electro['beta'] = `self._beta`
            
        for kf in self._kf:

            unit_fctr = (math.pow(_length[_ulen], -ldim) *
                         math.pow(_moles[_umol], -mdim) / _time[_utime])
            
            if type(kf) == types.InstanceType:
                k = kf
            else:
                k = Arrhenius(A = kf[0], n = kf[1], E = kf[2])
            k.build(kfnode, unit_fctr, gas_species = self._igspecies,
                    name = nm, rxn_phase = self._rxnphase)

            # set values for low-pressure rate coeff if falloff rxn
            mdim += 1
            ldim -= 3
            nm = 'k0'

        str = str[:-1]
        r.addChild('reactants',str)
        str = ''
        for s in self._p.keys():
            ns = self._p[s]
            str += s+':'+`ns`+' '
        str = str[:-1]
        r.addChild('products',str)
        return r

#-------------------


class three_body_reaction(reaction):
    def __init__(self,
                 equation = '',
                 kf = None,
                 efficiencies = '',
                 id = '',
                 options = []
                 ):

        reaction.__init__(self, equation, kf, id, '', options)
        self._type = 'threeBody'
        self._effm = 1.0
        self._eff = efficiencies

        # clean up reactant and product lists
        for r in self._r.keys():
            if r == 'M' or r == 'm':
                del self._r[r]
        for p in self._p.keys():
            if p == 'M' or p == 'm':
                del self._p[p]                

        
    def build(self, p):
        r = reaction.build(self, p)
        if r == 0: return
        kfnode = r.child('rateCoeff')
        
        if self._eff:
            eff = kfnode.addChild('efficiencies',self._eff)
            eff['default'] = `self._effm`
        

#---------------


class falloff_reaction(reaction):
    def __init__(self,
                 equation = '',
                 kf0 = None,
                 kf = None, 
                 efficiencies = '',
                 falloff = None,
                 id = '',
                 options = []
                 ):

        kf2 = (kf, kf0)        
        reaction.__init__(self, equation, kf2, id, '', options)
        self._type = 'falloff'
        # use a Lindemann falloff function by default
        self._falloff = falloff
        if self._falloff == None:
            self._falloff = Lindemann()
        
        self._effm = 1.0
        self._eff = efficiencies        

        # clean up reactant and product lists
        del self._r['(+']
        del self._p['(+']
        if self._r.has_key('M)'):
            del self._r['M)']
            del self._p['M)']
        if self._r.has_key('m)'):
            del self._r['m)']
            del self._p['m)']
        else:
            for r in self._r.keys():
                if r[-1] == ')' and r.find('(') < 0:
                    if self._eff:
                        raise CTI_Error('(+ '+mspecies+') and '+self._eff+' cannot both be specified')
                    self._eff = r[-1]+':1.0'
                    self._effm = 0.0
                    
                    del self._r[r]
                    del self._p[r]

        
    def build(self, p):
        r = reaction.build(self, p)
        if r == 0: return
        kfnode = r.child('rateCoeff')
        
        if self._eff and self._effm >= 0.0:
            eff = kfnode.addChild('efficiencies',self._eff)
            eff['default'] = `self._effm`
        
        if self._falloff:
            self._falloff.build(kfnode)


class surface_reaction(reaction):

    def __init__(self,
                 equation = '',
                 kf = None,
                 id = '',
                 order = '',                 
                 options = []):
        reaction.__init__(self, equation, kf, id, order, options)
        self._type = 'surface'


class edge_reaction(reaction):

    def __init__(self,
                 equation = '',
                 kf = None,
                 id = '',
                 order = '',
                 beta = 0.0,
                 options = []):
        reaction.__init__(self, equation, kf, id, order, options)
        self._type = 'edge'
        self._beta = beta


#--------------     


class state:
    def __init__(self,
                 temperature = None,
                 pressure = None,
                 mole_fractions = None,
                 mass_fractions = None,
                 density = None,
                 coverages = None):
        self._t = temperature
        self._p = pressure
        self._rho = density
        self._x = mole_fractions
        self._y = mass_fractions
        self._c = coverages
        
    def build(self, ph):
        st = ph.addChild('state')
        if self._t: addFloat(st, 'temperature', self._t, defunits = 'K')
        if self._p: addFloat(st, 'pressure', self._p, defunits = _upres)
        if self._rho: addFloat(st, 'density', self._rho, defunits = _umass+'/'+_ulen+'3')
        if self._x: st.addChild('moleFractions', self._x)
        if self._y: st.addChild('massFractions', self._y)
        if self._c: st.addChild('coverages', self._c)
    

class phase:
    """Base class for phases of matter."""
    
    def __init__(self,
                 name = '',
                 dim = 3,
                 elements = '',
                 species = '',
                 reactions = 'none',
                 initial_state = None,
                 options = []):
        
        self._name = name
        self._dim = dim
        self._el = elements
        self._sp = []
        self._rx = []
        
        if type(options) == types.StringType:
            self._options = [options]
        else:
            self._options = options        

        self.debug = 0
        if 'debug' in options:
            self.debug = 1

        #--------------------------------
        #        process species
        #--------------------------------
        
        # if a single string is entered, make it a list
        if type(species) == types.StringType:
            self._species = [species]
        else:
            self._species = species

        self._skip = 0

        # dictionary of species names 
        self._spmap = {}

        # for each species string, check whether or not the species
        # are imported or defined locally. If imported, the string
        # contains a colon (:) 
        for sp in self._species:
            icolon = sp.find(':')
            if icolon > 0:
                #datasrc, spnames = sp.split(':')
                datasrc = sp[:icolon].strip()
                spnames = sp[icolon+1:]
                self._sp.append((datasrc+'.xml', spnames))
            else:
                spnames = sp
                self._sp.append(('', spnames))                

            # strip the commas, and make the list of species names
            # 10/31/03: commented out the next line, so that species names may contain commas
            #sptoks = spnames.replace(',',' ').split()
            sptoks = spnames.split()
            
            for s in sptoks:
                # check for stray commas
                if s <> ',':
                    if s[0] == ',': s = s[1:]
                    if s[-1] == ',': s = s[:-1]

                    if self._spmap.has_key(s):
                        raise CTI_Error('Multiply-declared species '+s+' in phase '+self._name)
                    self._spmap[s] = self._dim
            
        self._rxns = reactions

        # check that species have been declared
        if len(self._spmap) == 0:
            raise CTI_Error('No species declared for phase '+self._name)
        
        # and that only one species is declared if it is a pure phase
        if self.is_pure() and len(self._spmap) > 1:
            raise CTI_Error('Stoichiometric phases must  declare exactly one species, \n'+
                               'but phase '+self._name+' declares '+`len(self._spmap)`+'.')

        self._initial = initial_state

        # add this phase to the global phase list
        global _phases
        _phases.append(self)
        

    def is_ideal_gas(self):
        """True if the entry represents an ideal gas."""
        return 0
    
    def is_pure(self):
        return 0

    def has_species(self, s):
        """Return 1 is a species with name 's' belongs to the phase,
        or 0 otherwise."""
        if self._spmap.has_key(s): return 1
        return 0

    def conc_dim(self):
        """Concentration dimensions. Used in computing the units for reaction
        rate coefficients."""
        return (1, -self._dim) 

    
    def buildrxns(self, p):
        
        if type(self._rxns) == types.StringType:
            self._rxns = [self._rxns]
            
        # for each reaction string, check whether or not the reactions
        # are imported or defined locally. If imported, the string
        # contains a colon (:) 
        for r in self._rxns:
            icolon = r.find(':')
            if icolon > 0:
                #datasrc, rnum = r.split(':')
                datasrc = r[:icolon].strip()
                rnum = r[icolon+1:]
                self._rx.append((datasrc+'.xml', rnum))
            else:
                rnum = r
                self._rx.append(('', rnum))                

        for r in self._rx:
            datasrc = r[0]
            ra = p.addChild('reactionArray')
            ra['datasrc'] = datasrc+'#reaction_data'
            if 'skip_undeclared_species' in self._options:
                rk = ra.addChild('skip')
                rk['species'] = 'undeclared'
            
            rtoks = r[1].split()
            if rtoks[0] <> 'all':
                i = ra.addChild('include')
                #i['prefix'] = 'reaction_'
                i['min'] = rtoks[0]
                if len(rtoks) > 2 and (rtoks[1] == 'to' or rtoks[1] == '-'):
                    i['max'] = rtoks[2]
                else:
                    i['max'] = rtoks[0]                
                    
            
    def build(self, p):
        p.addComment('    phase '+self._name+'     ')
        ph = p.addChild('phase')
        ph['id'] = self._name
        ph['dim'] = `self._dim`

        # ------- error tests -------
        #err = ph.addChild('validation')
        #err.addChild('duplicateReactions','halt')
        #err.addChild('thermo','warn')
        
        e = ph.addChild('elementArray',self._el)
        e['datasrc'] = 'elements.xml'
        for s in self._sp:
            datasrc, names = s
            sa = ph.addChild('speciesArray',names)
            sa['datasrc'] = datasrc+'#species_data'

            if 'skip_undeclared_elements' in self._options:
                sk = sa.addChild('skip')
                sk['element'] = 'undeclared'
            
        if self._rxns <> 'none':
            self.buildrxns(ph)
            
        #self._eos.build(ph)
        if self._initial:
            self._initial.build(ph)
        return ph


class ideal_gas(phase):
    """An ideal gas mixture."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 reactions = 'none',
                 kinetics = 'GasKinetics',
                 transport = 'None',
                 initial_state = None,
                 options = []):

        phase.__init__(self, name, 3, elements, species, reactions,
                       initial_state, options)
        self._pure = 0
        self._kin = kinetics
        self._tr = transport
        if self.debug:
            print 'Read ideal_gas entry '+self._name
            try:
                print 'in file '+__name__
            except:
                pass
                

        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'IdealGas'
        k = ph.addChild("kinetics")
        k['model'] = self._kin
        t = ph.addChild('transport')
        t['model'] = self._tr

    def is_ideal_gas(self):
        return 1


class stoichiometric_solid(phase):
    """A solid compound or pure element.Stoichiometric solid phases
    contain exactly one species, which always has unit activity. The
    solid is assumed to have constant density. Therefore the rates of
    reactions involving these phases do not contain any concentration
    terms for the (one) species in the phase, since the concentration
    is always the same.  """
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None,
                 options = []):
        
        phase.__init__(self, name, 3, elements, species, 'none',
                       initial_state, options)
        self._dens = density
        self._pure = 1
        if self._dens < 0.0:
            raise CTI_Error('density must be specified.')
        self._tr = transport

    def conc_dim(self):
        return (0,0)
        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'StoichSubstance'
        addFloat(e, 'density', self._dens, defunits = _umass+'/'+_ulen+'3')
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'        


class stoichiometric_liquid(stoichiometric_solid):
    """A stoichiometric liquid. Currently, there is no distinction
    between stoichiometric liquids and solids."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None,
                 options = []):
        
        stoichiometric_solid.__init__(self, name, 3, elements,
                                      species, 'none',
                                      initial_state, options)
        self._dens = density
        self._pure = 1
        if self._dens < 0.0:
            raise CTI_Error('density must be specified.')
        self._tr = transport
        
class pure_solid(stoichiometric_solid):
    """Deprecated. Use stoichiometric_solid"""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None,
                 options = []):
        
        stoichiometric_solid.__init__(self, name, 3, elements,
                                      species, 'none',
                                      initial_state, options)
        self._dens = density
        self._pure = 1
        if self._dens < 0.0:
            raise CTI_Error('density must be specified.')
        self._tr = transport
        print 'WARNING: entry type pure_solid is deprecated.'
        print 'Use stoichiometric_solid instead.'


class metal(phase):
    """A metal."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None,
                 options = []):
        
        phase.__init__(self, name, 3, elements, species, 'none',
                       initial_state, options)
        self._dens = density
        self._pure = 0
        self._tr = transport

    def conc_dim(self):
        return (0,0)
        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'Metal'
        addFloat(e, 'density', self._dens, defunits = _umass+'/'+_ulen+'3')
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'        


class incompressible_solid(phase):
    """An incompressible solid."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 density = -1.0,
                 transport = 'None',
                 initial_state = None,
                 options = []):
        
        phase.__init__(self, name, 3, elements, species, 'none',
                       initial_state, options)
        self._dens = density
        self._pure = 0
        if self._dens < 0.0:
            raise CTI_Error('density must be specified.')
        self._tr = transport

    def conc_dim(self):
        return (1,-3)
        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'Incompressible'
        addFloat(e, 'density', self._dens, defunits = _umass+'/'+_ulen+'3')
        if self._tr:
            t = ph.addChild('transport')
            t['model'] = self._tr
        k = ph.addChild("kinetics")
        k['model'] = 'none'        


class liquid_vapor(phase):
    """A fluid with a complete liquid/vapor equation of state.
    This entry type selects one of a set of predefined fluids with
    built-in liquid/vapor equations of state. The substance_flag
    parameter selects the fluid. See purefluids.py for the usage
    of this entry type."""
    
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 substance_flag = 0,
                 initial_state = None,
                 options = []):
        
        phase.__init__(self, name, 3, elements, species, 'none',
                       initial_state, options)
        self._subflag = substance_flag
        self._pure = 1


    def conc_dim(self):
        return (0,0)
        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'PureFluid'
        e['fluid_type'] = `self._subflag`
        k = ph.addChild("kinetics")
        k['model'] = 'none'        

class ideal_interface(phase):
    """An ideal interface."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 reactions = 'none',
                 site_density = 0.0,
                 phases = [],
                 kinetics = 'Interface',
                 transport = 'None',
                 initial_state = None,
                 options = []):

        self._type = 'surface'
        phase.__init__(self, name, 2, elements, species, reactions,
                       initial_state, options)
        self._pure = 0
        self._kin = kinetics
        self._tr = transport
        self._phases = phases
        self._sitedens = site_density
        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'Surface'
        addFloat(e, 'site_density', self._sitedens, defunits = _umol+'/'+_ulen+'2')
        k = ph.addChild("kinetics")
        k['model'] = self._kin
        t = ph.addChild('transport')
        t['model'] = self._tr
        p = ph.addChild('phaseArray',self._phases)


    def conc_dim(self):
        return (1, -2)


class edge(phase):
    """A 1D boundary between two surface phases."""
    def __init__(self,
                 name = '',
                 elements = '',
                 species = '',
                 reactions = 'none',
                 site_density = 0.0,
                 phases = [],
                 kinetics = 'Edge',
                 transport = 'None',
                 initial_state = None,
                 options = []):

        self._type = 'edge'
        phase.__init__(self, name, 1, elements, species, reactions,
                       initial_state, options)
        self._pure = 0
        self._kin = kinetics
        self._tr = transport
        self._phases = phases
        self._sitedens = site_density
        
    def build(self, p):
        ph = phase.build(self, p)
        e = ph.addChild("thermo")
        e['model'] = 'Edge'
        addFloat(e, 'site_density', self._sitedens, defunits = _umol+'/'+_ulen)
        k = ph.addChild("kinetics")
        k['model'] = self._kin
        t = ph.addChild('transport')
        t['model'] = self._tr
        p = ph.addChild('phaseArray',self._phases)


    def conc_dim(self):
        return (1, -1)    

        
#-------------------------------------------------------------------

# falloff parameterizations

class Troe:
    
    def __init__(self, A = 0.0, T3 = 0.0, T1 = 0.0, T2 = -999.9):
        if T2 <> -999.9:
            self._c = (A, T3, T1, T2)
        else:
            self._c = (A, T3, T1)
            
    def build(self, p):
        s = ''
        for num in self._c:
            s += '%g ' % num
        f = p.addChild('falloff', s)
        f['type'] = 'Troe'


class SRI:
    def __init__(self, A = 0.0, B = 0.0, C = 0.0, D = -999.9, E=-999.9):
        if D <> -999.9 and E <> -999.9:
            self._c = (A, B, C, D, E)
        else:
            self._c = (A, B, C)
        
    def build(self, p):
        s = ''
        for num in self._c:
            s += '%g ' % num
        f = p.addChild('falloff', s)
        f['type'] = 'SRI'


class Lindemann:
    def __init__(self):
        pass
    def build(self, p):
        f = p.addChild('falloff')
        f['type'] = 'Lindemann'

#--------------------------------------------------------------------

## class gas_transport:
##     def __init__(self, geom = 'nonlin',
##                  welldepth = 0.0,
##                  diam = 0.0,
##                  dipole = 0.0,
##                  polar = 0.0,
##                  rot_relax = 0.0):
##         self._sp = species
##         self._geom = geom
##         self._params = (welldepth, diam, dipole, polar, rotrelax)

##         #global _trdata
##         #_trdata[species] = self

##     def build(self, s):
##         tr = s.addChild('transport')
##         g = tr.addChild('string','linear')
##         g['title'] = 'geometry'
##         tr.addChild('LJ_welldepth',`self._params[0]`)
##         tr.addChild('LJ_diameter',`self._params[1]`)
##         tr.addChild('dipoleMoment',`self._params[2]`)
##         tr.addChild('polarizability',`self._params[3]`)
##         tr.addChild('rotRelax',`self._params[4]`)        


#get_atomic_wts()
validate()


## if __name__ == "__main__":
##     from Cantera import *
##     import sys, os, os.path
##     file = sys.argv[1]
##     base = os.path.basename(file)
##     root, ext = os.path.splitext(base)
##     dataset(root)
##     execfile(file)
##     write()
    

##########################################
#
# $Author$
# $Revision$
# $Date$
# $Log$
# Revision 1.5  2004-09-20 10:25:02  dggoodwin
# *** empty log message ***
#
# Revision 1.4  2004/07/14 11:24:13  dggoodwin
# *** empty log message ***
#
# Revision 1.3  2004/06/09 01:02:31  dggoodwin
# *** empty log message ***
#
# Revision 1.2  2004/06/04 06:05:31  dggoodwin
# *** empty log message ***
#
# Revision 1.1  2004/06/02 04:39:09  dggoodwin
# moved ctml_writer.py out of Cantera package
#
# Revision 1.34  2004/05/30 04:02:55  dggoodwin
# converted to pure python
#
# Revision 1.33  2004/04/23 19:03:21  dggoodwin
# *** empty log message ***
#
# Revision 1.32  2004/04/23 16:35:32  dggoodwin
# *** empty log message ***
#
# Revision 1.31  2004/03/12 05:59:59  dggoodwin
# *** empty log message ***
#
# Revision 1.30  2004/02/08 13:25:21  dggoodwin
# *** empty log message ***
#
# Revision 1.29  2004/02/08 13:22:31  dggoodwin
# *** empty log message ***
#
# Revision 1.28  2004/02/08 13:09:10  dggoodwin
# *** empty log message ***
#
# Revision 1.27  2004/02/03 16:42:54  dggoodwin
# *** empty log message ***
#
# Revision 1.26  2004/02/03 03:31:06  dggoodwin
# *** empty log message ***
#
# Revision 1.25  2003/11/24 16:39:33  dggoodwin
# -
#
# Revision 1.24  2003/11/13 12:29:45  dggoodwin
# *** empty log message ***
#
# Revision 1.23  2003/11/12 18:58:15  dggoodwin
# *** empty log message ***
#
# Revision 1.22  2003/11/01 04:48:20  dggoodwin
# added capability to have species names with embedded commas
#
# Revision 1.21  2003/10/14 06:48:07  dggoodwin
# *** empty log message ***
#
# Revision 1.20  2003/09/22 13:14:34  dggoodwin
# *** empty log message ***
#
# Revision 1.19  2003/08/26 03:39:02  dggoodwin
# *** empty log message ***
#
# Revision 1.18  2003/08/21 14:29:53  dggoodwin
# *** empty log message ***
#
# Revision 1.17  2003/08/20 15:35:32  dggoodwin
# *** empty log message ***
#
# Revision 1.16  2003/08/19 22:02:01  hkmoffa
# Fixed an error in an argument list
#
# Revision 1.15  2003/08/18 05:05:02  dggoodwin
# added support for specified reaction order, sticking coefficients,
# coverage dependence of rate coefficients; fixed error where site_density
# was not being converted to SI.
#
# Revision 1.14  2003/08/16 20:17:21  dggoodwin
# changed handling of reaction pre-exponential units to write converted
# value to CTML, rather than pass original value with a units string
#
#
###########################################
