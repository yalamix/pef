from typing import List
from copy import copy
from .graph import *
from .utils import *
from sympy import *
import json

BEAM_HEIGHT = 0.1
BEAM_COLOR = 'blue'
DEFAULT_BEAM = {
    "beam_size": 5,
    "variables": {
        "L": 5,
        "A": 1,
        "E": 200 * 1e9
    },
    "points": {
        "A": 0,
        "B": 'L'
    },
    "links": [],
    "shear_forces": [],
    "normal_forces": [],
    "twisting_moments": [],
    "bending_moments": []
}

class BeamProblem:
    def __init__(self, definition: str):
        """
        Creates a BeamProblem from a json with its definitions.
        
        Parameters:
            definition         : str of a json.
        
        Returns:
            The BeamProblem object.
        """           
        b = json.loads(definition)
        self.beam_size = b['beam_size']
        self.variables = b['variables']
        self.points = b['points']
        self.links = b['links']
        self.shear_forces = b['shear_forces']
        self.normal_forces = b['normal_forces']
        self.twisting_moments = b['twisting_moments']
        self.bending_moments = b['bending_moments']
        self.protected_symbols = [
            'N_x',
            'M_x',
            'M_z',
            'V_y',
            'u',
            'q',
            'p',
            'v',
            'phi',
            'theta_Z',
            'tau',
            'E',
            'G',
            'A',
            'J_p',
            'I_zz'
        ]
        self.symbols = {}
        self._update_symbols()
    
    def to_dict(self) -> dict:
        """
        Returns a dictionary with the problem definitions.
        
        Returns:
            b: The problem dictionary.
        """              
        b = {}
        b['beam_size'] = self.beam_size
        b['variables'] = self.variables
        b['points'] = self.points
        b['links'] = self.links
        b['shear_forces'] = self.shear_forces
        b['normal_forces'] = self.normal_forces
        b['twisting_moments'] = self.twisting_moments
        b['bending_moments'] = self.bending_moments
        return b
    
    def _update_variables(self, expr: str, default: float = 1) -> str:
        """
        Updates the variable list.   

        Parameters:
            expr        : Expression containing new variables.
            default        : Default value for new variables.
        """          
        new_symbols = parse_and_update_symbols(expr, list(self.variables.keys()))
        for s in new_symbols:
            if s not in self.variables and s not in self.protected_symbols:
                self.variables[s] = default
            if s not in self.variables and s in self.protected_symbols:
                sn = s + '_1'
                self.variables[sn] = default
                expr = expr.replace(s,sn)
        self._update_symbols()
        return str(sympify(expr))
    
    def _update_symbols(self):
        for name in list(self.variables.keys()) + self.protected_symbols:
            if name not in self.symbols:
                self.symbols[name] = Symbol(name)
    
    def _remove_protected_symbols(self, v: dict) -> dict:
        b = {}
        for item in v:
            if str(item) not in self.protected_symbols:
                b[item] = v[item]
        return b

    def _get_symbol_value(self) -> dict:
        """
        Gets Symbol -> value pairs.

        Returns:
            Symbol dictionary.
        """
        b = {}
        sym, f = self._get_functions_and_symbols()
        for s in sym:
            if s in self.variables:
                b[sym[s]] = self.variables[s]
        return b

    def _get_functions_and_symbols(self):
        """
        Gets all functions and symbols
        """
        s = self.symbols
        f = {
            name: Function(s[name]) for name in [
                'N_x',
                'M_x',
                'M_z',
                'V_y',
                'u',
                'q',
                'p',
                'v',
                'phi',
                'theta_Z'
            ]
        }
        return s, f        
        
    def _s(self, expr: str):
        return sympify(expr, locals=self.symbols)

    def _ev(self, expr: str) -> float:
        """
        Evaluates a Sympy expression.

        Returns:
            Resulting value.
        """ 
        if isinstance(expr, float) or isinstance(expr, int):
            return expr
        return float(self._s(expr).evalf(subs=self.variables))
    
    def _get_and_add_point(self, pos: str) -> str:
        """
        Adds a new point of interest and returns it. If it already exists, return it.

        Returns:
            current: Newly added point.        
        """
        last_point = 'A'
        for point in self.points:
            if self._ev(self.points[point]) == self._ev(pos):
                return point
            last_point = point
            if self._ev(self.points[point]) > self._ev(pos):
                break
        current = last_point            
        old_points = copy(self.points)
        print()
        for point in old_points:
            if ord(point) == ord(last_point):
                self.points[last_point] = pos
                self.points[chr(ord(point) + 1)] = old_points[point]
                last_point = chr(ord(point) + 1)
                pos = self.points[last_point]
        return current
    
    def _enforce_reaction_consistency(self) -> None:
        """
        Enforces reaction naming convention.
        """
        for point in self.points:
            for force in self.shear_forces:
                if 'R_y_' in force['value'] and force['n'] < 0:
                    if self._ev(force['start']) == self._ev(self.points[point]):
                        force['value'] = f'R_y_{point}'

    def _assert_link_consistency(self, link_type: str, position: float) -> bool:
        """
        Makes sure that two links don't occupy the same position.

        Parameters:
            link_type        : Name of the link.        
            position         : Position of the link.        
        
        Returns:
            True if the link can be added, False otherwise.
        """
        for link in self.links:
            if link_type in link:
                try:
                    assert link[link_type] != position
                except:
                    return False
        return True
    
    def _enforce_correct_rescaling(self, old: float, new: float) -> None:
        """
        Rescales positions in case the user changes the beam length.
        """         
        def rescale(item, o, n):
            if isinstance(item, float) or isinstance(item, int):
                return (new/old) * item
            else:
                return f'({str(item)}) * ({new/old})'
        # loads
        for f in self.shear_forces:
            if self._ev(f['start']) > new:
                f['start'] = rescale(f['start'], old, new)
            if self._ev(f['stop']) > new:
                f['stop'] = rescale(f['stop'], old, new)
        for f in self.normal_forces:
            if self._ev(f['start']) > new:
                f['start'] = rescale(f['start'], old, new)
            if self._ev(f['stop']) > new:
                f['stop'] = rescale(f['stop'], old, new)      
        for f in self.twisting_moments:
            if self._ev(f['start']) > new:
                f['start'] = rescale(f['start'], old, new)
            if self._ev(f['stop']) > new:
                f['stop'] = rescale(f['stop'], old, new)     
        for f in self.bending_moments:
            if self._ev(f['start']) > new:
                f['start'] = rescale(f['start'], old, new)
            if self._ev(f['stop']) > new:
                f['stop'] = rescale(f['stop'], old, new)
        # links
        for link in self.links:
            if self._ev(link[list(link.keys())[0]]) > new:
                link[list(link.keys())[0]] = rescale(link[list(link.keys())[0]], old, new)
    
    def remove_item(self, item_type: str, item_index: int):
        """
        Removes an item (link, variable, force, moment) from the problem.   

        Parameters:
            item_type        : Type of item to remove.        
            item_index         : Position of item in its item type list.
        """ 
        if item_type == 'link':
            self.links.pop(item_index)
        if item_type == 'variable':
            self.variables.pop(item_index)
        if item_type == 'shear_force':
            self.shear_forces.pop(item_index)
        if item_type == 'normal_force':
            self.normal_forces.pop(item_index)
        if item_type == 'twisting_moment':
            self.twisting_moments.pop(item_index)
        if item_type == 'bending_moment':
            self.bending_moments.pop(item_index)         

    def add_variable(self, variable: str, value: float):
        """
        Adds a new variable. If already exists, updates value.

        Parameters:
            variable        : Variable name/symbol.        
            value         : Variable value.
        """ 
        self.variables[variable] = value
        if variable == 'L':
            self._enforce_correct_rescaling(self.beam_size, value)
            self.beam_size = value
          
    def add_link(self, link_type: str, position: str) -> bool:
        """
        Adds a new link.

        Parameters:
            link_type        : Name of the link.        
            position         : Position of the link.             
        
        Returns:
            True if the link was added, False otherwise.
        """
        # check position
        try:
            position = float(position)
        except:
            if position not in self.variables:
                position = self._update_variables(position)

        if self._assert_link_consistency(link_type, position):
            # enforces cantilevers being placed at the boundaries
            if link_type == 'cantilever':
                if self._ev(position) < self.beam_size/2:
                    position = 0
                else:
                    position = 'L'
            # enforce that links be placed inside problem domain
            if self._ev(position) < 0:
                position = 0
            if self._ev(position) > self.beam_size:
                position = 'L'
            # creates point of interest
            if link_type in ('fixed_support', 'mobile_support') and (0 < self._ev(position) < self.beam_size):
                link_point = self._get_and_add_point(position)
                self.add_shear_force(f'R_y_{link_point}', position, position, -1, True, True)
                self._enforce_reaction_consistency()
            else:
                link_point = self._get_and_add_point(position)
                self._enforce_reaction_consistency()
            self.links.append({link_type: position})
            return True
        return False

    def _create_load(self, value: str, start: str, stop: str, n: int, pos: bool = True, reaction: bool = False) -> dict:
        """
        Creates a load (force or moment).   

        Parameters:
            value        : Value of the force.        
            start         : Start position of the domain of the force.
            stop         : End position of the domain of the force.
            n         : Exponent of the polynomial representing the force.
            pos         : True if the force is positive, False otherwise.
        """         
        if value not in self.variables:
            try:
                float(value)
            except:
                if reaction:
                    if reaction not in self.protected_symbols:
                        self.protected_symbols.append(value)
                else:
                    value = self._update_variables(value)
        if start not in self.variables:
            try:
                float(start)
            except:     
                start = self._update_variables(start)
        if stop not in self.variables:
            try:
                float(stop)
            except:     
                stop = self._update_variables(stop, 2)

        # if n >= 0 assert stop > start
        if n >= 0:
            if self._ev(stop) < self._ev(start):
                self.stop = 'L'
        # constrain load to problem domain
        if self._ev(start) < 0:
            start = 0
        if self._ev(stop) > self.beam_size:
            stop = 'L'

        # add points of interest
        self._get_and_add_point(start)
        if self._ev(start) != self._ev(stop):
            self._get_and_add_point(stop)
        self._enforce_reaction_consistency()

        return {
                "value": value,
                "start": start,
                "stop": stop,
                "n": n,
                "pos": pos,
                "r": reaction
            }          

    def add_shear_force(self, value: str, start: str, stop: str, n: int, pos: bool = True, reaction: bool = False):
        """
        Adds a concentrated or distributed shear force.   

        Parameters:
            value        : Value of the force.        
            start         : Start position of the domain of the force.
            stop         : End position of the domain of the force.
            n         : Exponent of the polynomial representing the force.
            pos         : True if the force is positive, False otherwise.
        """ 
        self.shear_forces.append(self._create_load(value, start, stop, n, pos, reaction))  

    def add_normal_force(self, value: str, start: str, stop: str, n: int, pos: bool = True):
        """
        Adds a concentrated or distributed normal force.   

        Parameters:
            value        : Value of the force.        
            start         : Start position of the domain of the force.
            stop         : End position of the domain of the force.
            n         : Exponent of the polynomial representing the force.
            pos         : True if the force is positive, False otherwise.
        """       
        self.normal_forces.append(self._create_load(value, start, stop, n, pos))  

    def add_twisting_moment(self, value: str, start: str, stop: str, n: int, pos: bool = True):
        """
        Adds a concentrated or distributed twisting moment.   

        Parameters:
            value        : Value of the moment.        
            start         : Start position of the domain of the moment.
            stop         : End position of the domain of the moment.
            n         : Exponent of the polynomial representing the moment.
            pos         : True if the moment is positive, False otherwise.
        """        
        self.twisting_moments.append(self._create_load(value, start, stop, n, pos))  

    def add_bending_moment(self, value: str, start: str, stop: str, n: int, pos: bool = True):
        """
        Adds a concentrated or distributed bending moment.   

        Parameters:
            value        : Value of the moment.        
            start         : Start position of the domain of the moment.
            stop         : End position of the domain of the moment.
            n         : Exponent of the polynomial representing the moment.
            pos         : True if the moment is positive, False otherwise.
        """        
        self.bending_moments.append(self._create_load(value, start, stop, n, pos))
    
    def _get_load_equation(self, force: dict):
        x = symbols('x')
        try:
            p = float(force['value'])
        except:
            p = self._s(force['value'])
        if not force['pos']:
            p = -p
        if self._ev(force['stop']) < self.beam_size and force['n'] >= 0:
            if force['n'] <= 0:
                eq = p * SingularityFunction(x, self._s(force['start']), force['n']) - p * SingularityFunction(x, self._s(force['stop']), force['n'])
            else:
                eq = Rational(p,2,gcd=1) * SingularityFunction(x, self._s(force['start']), force['n']) - Rational(p,2,gcd=1) * SingularityFunction(x, self._s(force['stop']), force['n'])
        else:
            if force['n'] <= 0:
                eq = p * SingularityFunction(x, self._s(force['start']), force['n'])
            else:
                eq = Rational(p,2,gcd=1) * SingularityFunction(x, self._s(force['start']), force['n'])
        return eq
    
    def _get_load_expression(self, force: dict) -> str:
        """
        Generates the latex string of the singularity function representing the load.   

        Parameters:
            force        : Load dictionary.        

        Returns:
            out: Latex string of the load expression.
        """          
        return latex_with_threshold(self._get_load_equation(force))
    
    def _load_expressions(self, loads: List[dict]) -> List[str]:
        """
        Returns a list of the latex form of the loads.
        
        Returns:
            b: The list of latex strings.
        """        
        b = []
        for load in loads:
            if load['n'] < 0 and (self._ev(load['start']) == 0 or self._ev(load['start'] == self.beam_size)):
                continue
            else:
                b.append(self._get_load_expression(load))
        return b

    def get_normal_forces(self) -> List[str]:
        """
        Returns a list of the latex form of the normal forces.
        
        Returns:
            b: The list of latex strings.
        """
        return self._load_expressions(self.normal_forces)

    def get_shear_forces(self) -> List[str]:
        """
        Returns a list of the latex form of the shear forces.
        
        Returns:
            b: The list of latex strings.
        """
        return self._load_expressions(self.shear_forces) 
 
    def get_twisting_moments(self) -> List[str]:
        """
        Returns a list of the latex form of the twisting moments.
        
        Returns:
            b: The list of latex strings.
        """
        return self._load_expressions(self.twisting_moments)      

    def get_bending_moments(self) -> List[str]:
        """
        Returns a list of the latex form of the bending moments.
        
        Returns:
            b: The list of latex strings.
        """
        return self._load_expressions(self.bending_moments)

    def get_boundary_conditions(self) -> List[str]:
        """
        Returns a list of the latex form of the boundary conditions.
        
        Returns:
            b: The list of latex strings.
        """ 
        bc, mp = self.calculate_boundary_conditions()
        return [latex_with_threshold(eq) for eq in bc]

    def _check_boundaries(self, left=True) -> bool:
        """
        """
        t = 0 if left else self.beam_size
        for load in self.normal_forces + self.shear_forces + self.bending_moments + self.twisting_moments:
            if self._ev(load['start']) == t:
                return True
        for link in self.links:
            if self._ev(link[list(link.keys())[0]]) == t:
                return True
        return False

    def calculate_boundary_conditions(self):
        """
        """ 
        conditions = []
        s, f = self._get_functions_and_symbols()
        mapped_points = {name: [] for name in [str(fun) for fun in f]}

        degrees_of_freedom = {
            'shear': 4,
            'normal': 2,
            'twisting': 2
        }

        # Normal
        for force in self.normal_forces:
            if force['r']:
                degrees_of_freedom['normal'] += 1
            if force['n'] < 0:
                if self._ev(force['start']) == 0:
                    if force['pos']:
                        conditions.append(Eq(f['N_x'](0),self._s(force['value'])))
                        mapped_points['N_x'].append(0)
                    else:
                        conditions.append(Eq(f['N_x'](0),-self._s(force['value'])))
                        mapped_points['N_x'].append(0)
                elif self._ev(force['start']) == self.beam_size:
                    if force['pos']:
                        conditions.append(Eq(f['N_x'](s['L']),self._s(force['value'])))
                        mapped_points['N_x'].append('L')
                    else:
                        conditions.append(Eq(f['N_x'](s['L']),-self._s(force['value'])))
                        mapped_points['N_x'].append('L')

        # Twisting
        for moment in self.twisting_moments:
            if moment['r']:
                degrees_of_freedom['twisting'] += 1            
            if moment['n'] < 0:
                if self._ev(moment['start']) == 0:
                    if moment['pos']:
                        conditions.append(Eq(f['M_x'](0),self._s(force['value'])))
                        mapped_points['M_x'].append(0)
                    else:
                        conditions.append(Eq(f['M_x'](0),-self._s(force['value'])))
                        mapped_points['M_x'].append(0)
                elif self._ev(moment['start']) == self.beam_size:
                    if moment['pos']:
                        conditions.append(Eq(f['M_x'](s['L']),self._s(force['value'])))
                        mapped_points['M_x'].append('L')
                    else:
                        conditions.append(Eq(f['M_x'](s['L']),-self._s(force['value'])))
                        mapped_points['M_x'].append('L')  

        # Shear
        for force in self.shear_forces:
            if force['r']:
                degrees_of_freedom['shear'] += 1            
            if force['n'] < 0:
                if self._ev(force['start']) == 0:
                    if force['pos']:
                        conditions.append(Eq(f['V_y'](0),self._s(force['value'])))
                        mapped_points['V_y'].append(0)
                    else:
                        conditions.append(Eq(f['V_y'](0),-self._s(force['value'])))
                        mapped_points['V_y'].append(0)
                elif self._ev(force['start']) == self.beam_size:
                    if force['pos']:
                        conditions.append(Eq(f['V_y'](s['L']),self._s(force['value'])))
                        mapped_points['V_y'].append('L')  
                    else:
                        conditions.append(Eq(f['V_y'](s['L']),-self._s(force['value'])))
                        mapped_points['V_y'].append('L')   

        # Bending
        for moment in self.bending_moments:
            if moment['r']:
                degrees_of_freedom['shear'] += 1               
            if moment['n'] < 0:
                if self._ev(moment['start']) == 0:
                    if moment['pos']:
                        conditions.append(Eq(f['M_z'](0),self._s(moment['value'])))
                        mapped_points['M_z'].append(0)
                    else:
                        conditions.append(Eq(f['M_z'](0),-self._s(moment['value'])))
                        mapped_points['M_z'].append(0)
                elif self._ev(moment['start']) == self.beam_size:
                    if moment['pos']:
                        conditions.append(Eq(f['M_z'](s['L']),self._s(moment['value'])))
                        mapped_points['M_z'].append('L')
                    else:
                        conditions.append(Eq(f['M_z'](s['L']),-self._s(moment['value'])))
                        mapped_points['M_z'].append('L')

        # Links
        for link in self.links:
            link_type = list(link.keys())[0]
            position = link[link_type]
            if link_type == 'cantilever':
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    # Shear and bending
                    conditions.append(Eq(f['V_y'](self._s(position)),0))
                    mapped_points['V_y'].append(position)                
                    conditions.append(Eq(f['M_z'](self._s(position)),0))
                    mapped_points['M_z'].append(position)
                    if len(mapped_points['M_z']) + len(mapped_points['V_y']) < degrees_of_freedom['shear']:
                        conditions.append(Eq(f['theta_Z'](self._s(position)),0))
                        mapped_points['theta_Z'].append(position)                  
                        conditions.append(Eq(f['v'](self._s(position)),0))
                        mapped_points['v'].append(position)      
                if len(self.normal_forces) > 0:
                    # Normal     
                    if len(mapped_points['N_x']) + len(mapped_points['u']) < degrees_of_freedom['normal']:
                        conditions.append(Eq(f['u'](self._s(position)),0))
                        mapped_points['u'].append(position)  
                if len(self.twisting_moments) > 0:
                    # Twisting                            
                    conditions.append(Eq(f['M_x'](self._s(position)),0))
                    mapped_points['M_x'].append(position)
                    if len(mapped_points['M_x']) + len(mapped_points['phi']) < degrees_of_freedom['twisting']:
                        conditions.append(Eq(f['phi'](self._s(position)),0))
                        mapped_points['phi'].append(position)                            
            if link_type in ('fixed_support', 'mobile_support'):
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    conditions.append(Eq(f['v'](self._s(position)),0))
                    mapped_points['v'].append(position)
            if link_type == 'hinge':
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    conditions.append(Eq(f['M_z'](self._s(position)),0))
                    mapped_points['M_z'].append(position)       
            if link_type == 'roller':
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    conditions.append(Eq(f['V_y'](self._s(position)),0))
                    mapped_points['V_y'].append(position)     

        if len(mapped_points['M_z']) + len(mapped_points['V_y']) + len(mapped_points['theta_Z']) + len(mapped_points['v']) < degrees_of_freedom['shear']:
            if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                if 0 not in mapped_points['V_y'] and not self._check_boundaries():
                    conditions.append(Eq(f['V_y'](0),0))
                    mapped_points['V_y'].append(0)
                if 'L' not in mapped_points['V_y'] and not self._check_boundaries(False):
                    conditions.append(Eq(f['V_y'](s['L']),0))
                    mapped_points['V_y'].append('L')                
                if 0 not in mapped_points['M_z'] and not self._check_boundaries():                
                    conditions.append(Eq(f['M_z'](0),0))
                    mapped_points['M_z'].append(0)
                if 'L' not in mapped_points['M_z'] and not self._check_boundaries(False):
                    conditions.append(Eq(f['M_z'](s['L']),0))
                    mapped_points['M_z'].append('L')              

        if len(mapped_points['N_x']) + len(mapped_points['u']) < degrees_of_freedom['normal']:
            if len(self.normal_forces) > 0:
                if 0 not in mapped_points['N_x'] and not self._check_boundaries():
                    conditions.append(Eq(f['N_x'](0),0))
                    mapped_points['N_x'].append(0)
                if 'L' not in mapped_points['N_x'] and not self._check_boundaries(False):
                    conditions.append(Eq(f['N_x'](s['L']),0))
                    mapped_points['N_x'].append('L')                  
        
        return conditions, mapped_points

    def solve(self):
        """
        """
        s, f = self._get_functions_and_symbols()
        bc, mp = self.calculate_boundary_conditions()
        x = symbols('x')
        vardict = self._get_symbol_value()
        solution_blocks = []
        # Normal
        normal_block = {'name': 'Força Normal'}
        reactions = []
        p = sympify(0)
        for force in self.normal_forces:
            if (force['n'] < 0 and (0 < self._ev(force['start']) < self.beam_size)) or force['n'] >= 0:
                p += self._get_load_equation(force)
            if force['r']:
                reactions.append(self._s(force['value']))
        normal_block['load_equation'] = latex_with_threshold(Eq(f['p'](x), p))
        c = []
        if len(mp['N_x']) > 0:
            if len(mp['u']) > 0:
                normal_block['differential_equation'] = latex_with_threshold(Eq(Eq(s['A']*s['E'] * Derivative(f['u'](x),(x,2)), -f['p'](x)), -p, evaluate=False))
                steps = []
                final_eqs = []
                plots = {}                
                c.append(Symbol('c_1'))
                p1 = integrate(-p, x) + c[0]
                steps.append(latex_with_threshold(Eq(Eq(s['A']*s['E'] * Derivative(f['u'](x),(x,1)),f['N_x'](x)), p1,evaluate=False)))
                final_eqs.append(Eq(f['N_x'](x), p1,evaluate=False))
                c.append(Symbol('c_2'))
                p2 = integrate(p1, x) + c[1]
                steps.append(latex_with_threshold(Eq(s['A']*s['E'] * f['u'](x),p2,evaluate=False)))
                final_eqs.append(Eq(s['A']*s['E'] * f['u'](x),p2,evaluate=False))
                normal_block['integration_steps'] = steps
                constant_det = {}
                constant_value = {}
                total_unknowns = len(c) + len(reactions)
                for cond in bc:
                    if cond.has(f['N_x']):
                        k = c[0]
                        p = p1
                        if str(k) in constant_det:
                            k = c[1]
                        constant_det[str(k)] = []
                        constant_det[str(k)].append(latex_with_threshold(Eq(cond, p, evaluate=False)))
                        if isinstance(cond.lhs.args[0], Symbol):
                            sb = Dummy(str(cond.lhs.args[0]), real=True, positive=True)
                            vardict[sb] = vardict[self.symbols[str(cond.lhs.args[0])]]
                            cond.xreplace({cond.lhs.args[0]: sb})
                        else:
                            sb = float(cond.lhs.args[0])
                        sol = solve(Eq(cond.rhs, p),k)[0]
                        print(sol, sb)
                        try:
                            sol = simplify(eval_all_singularities(Eq(k, sol), x, sb))
                        except:
                            sol = simplify(Eq(k, sol.subs(x, sb)))
                        print(sol)
                        constant_det[str(k)].append(latex_with_threshold(sol))
                        if latex_with_threshold(sol.evalf(subs=vardict)) != latex_with_threshold(sol):
                            constant_det[str(k)].append(latex_with_threshold(sol.evalf(subs=vardict)))
                        # sanity check
                        print("Before sub:", sol.rhs, "free symbols:", sol.rhs.free_symbols)
                        sym = list(sol.free_symbols)
                        if len(sym) > 0:
                            print("After sub:", sol.rhs.subs(sym[0], 5).evalf())
                        #
                        print("Before evalf:", sol.rhs, "free symbols:", sol.rhs.free_symbols)
                        print("Subs dictionary:", vardict)
                        print("After evalf:", sol.rhs.evalf(subs=vardict),'\n')   
                        constant_value[str(k)] = sol.rhs.evalf(subs=vardict)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                    if cond.has(f['u']):
                        k = c[1]
                        p = (1/(s['A']*s['E'])) * p2
                        constant_det[str(k)] = []
                        constant_det[str(k)].append(latex_with_threshold(Eq(cond, p, evaluate=False)))
                        if isinstance(cond.lhs.args[0], Symbol):
                            sb = Dummy(str(cond.lhs.args[0]), real=True, positive=True)
                            vardict[sb] = vardict[self.symbols[str(cond.lhs.args[0])]]
                            cond.xreplace({cond.lhs.args[0]: sb})
                        else:
                            sb = float(cond.lhs.args[0])
                        try:
                            sol = simplify(eval_all_singularities(Eq(k, solve(Eq(cond.rhs, p),k)[0]), x, sb))
                        except:
                            sol = simplify(Eq(k, solve(Eq(cond.rhs, p),k)[0].subs(x, sb)))
                        constant_det[str(k)].append(latex_with_threshold(sol))
                        if latex_with_threshold(sol.evalf(subs=vardict)) != latex_with_threshold(sol):
                            constant_det[str(k)].append(latex_with_threshold(sol.evalf(subs=vardict)))
                        # sanity check
                        print("Before sub:", sol.rhs, "free symbols:", sol.rhs.free_symbols)
                        sym = list(sol.free_symbols)
                        if len(sym) > 0:
                            print("After sub:", sol.rhs.subs(sym[0], 5).evalf())
                        #
                        print("Before evalf:", sol.rhs, "free symbols:", sol.rhs.free_symbols)
                        print("Subs dictionary:", vardict)
                        print("After evalf:", sol.rhs.evalf(subs=vardict),'\n')                       
                        constant_value[str(k)] = sol.rhs.evalf(subs=vardict)
                        if len(constant_det.keys()) == total_unknowns:
                            break                        
                normal_block['constants'] = constant_det
                for i in range(len(final_eqs)):
                    x_axis = np.linspace(0, self.variables['L'])
                    eq = final_eqs[i].subs(c[0],constant_value[str(c[0])]).subs(c[1],constant_value[str(c[1])]).evalf(subs=self._remove_protected_symbols(vardict))
                    f = lambdify(x, eq.rhs, modules=[mapping, 'numpy'])
                    plots[format_label(str(eq.lhs))] = fig_to_rotated_img(create_filled_line_figure(
                        x_axis,
                        f(x_axis),
                        format_label(str(eq.lhs))
                    ))
                    final_eqs[i] = latex_with_threshold(eq)
                normal_block['final_equations'] = final_eqs
                normal_block['plots'] = plots

            else:
                normal_block['differential_equation'] = latex_with_threshold(Eq(Eq(Derivative(f['N_x'](x)), -f['p'](x)), -p, evaluate=False))
                steps = []
                final_eqs = []
                plots = {}
                c.append(Symbol('c'))
                p = integrate(-p, x) + c[0]
                steps.append(latex_with_threshold(Eq(f['N_x'](x), p, evaluate=False)))
                final_eqs.append(Eq(f['N_x'](x), p, evaluate=False))
                normal_block['integration_steps'] = steps
                constant_det = {}
                constant_value = {}
                total_unknowns = len(c) + len(reactions)
                for cond in bc:
                    if cond.has(f['N_x']):
                        constant_det[str(c[0])] = []
                        constant_det[str(c[0])].append(latex_with_threshold(Eq(cond, p, evaluate=False)))
                        if isinstance(cond.lhs.args[0], Symbol):
                            sb = Dummy(str(cond.lhs.args[0]), real=True, positive=True)
                            vardict[sb] = vardict[self.symbols[str(cond.lhs.args[0])]]
                            cond.xreplace({cond.lhs.args[0]: sb})
                        else:
                            sb = float(cond.lhs.args[0])
                        try:
                            sol = simplify(eval_all_singularities(Eq(c[0], solve(Eq(cond.rhs, p),c[0])[0]), x, sb))
                        except:
                            sol = simplify(Eq(c[0], solve(Eq(cond.rhs, p),c[0])[0].subs(x, sb)))
                        constant_det[str(c[0])].append(latex_with_threshold(sol))
                        constant_det[str(c[0])].append(latex_with_threshold(sol.evalf(subs=vardict)))
                        constant_value[str(c[0])] = sol.rhs.evalf(subs=vardict)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                normal_block['constants'] = constant_det
                for i in range(len(final_eqs)):
                    x_axis = np.linspace(0, self.variables['L'])
                    eq = final_eqs[i].subs(c[0],constant_value[str(c[0])]).evalf(subs=vardict)
                    f = lambdify(x, eq.rhs, modules=[mapping, 'numpy'])
                    plots[format_label(str(eq.lhs))] = fig_to_rotated_img(create_filled_line_figure(
                        x_axis,
                        f(x_axis),
                        format_label(str(eq.lhs))
                    ))
                    final_eqs[i] = latex_with_threshold(eq)
                normal_block['final_equations'] = final_eqs
                normal_block['plots'] = plots
        
        solution_blocks.append(normal_block)

        return solution_blocks

    def graph(self) -> go.Figure:
        """
        Creates a Plotly figure of the Beam.
        
        Returns:
            fig (go.Figure): The Plotly figure.
        """         
        HEIGHT = (BEAM_HEIGHT * self.beam_size)
        # Create the initial figure of the beam
        self.fig = plot_rectangle(self.beam_size, HEIGHT, BEAM_COLOR)

        # Add the x and y axis arrows
        self.fig = add_axis_arrows(self.fig, HEIGHT)

        # Add points of interest
        total_points = len(list(self.points))
        for point in self.points:
            self.fig = add_label(self.fig, self._ev(self.points[point]), -2 * HEIGHT, point, font_size = 16 if total_points < 6 else 12)
            if ord(point) > ord('A') and total_points > 2:
                dist = str(simplify(self._s(f'({self.points[point]}) - ({self.points[chr(ord(point) - 1)]})'))).replace('*','')
                self.fig = add_hline_label(self.fig, -2 * HEIGHT, self._ev(self.points[chr(ord(point) - 1)]) + HEIGHT/8, self._ev(self.points[point]) - HEIGHT/8, format_label(dist))

        # Add the links
        for link in self.links:
            if 'cantilever' in link:
                self.fig = add_cantilever(self.fig, HEIGHT, self._ev(link['cantilever']))
            if 'hinge' in link:
                self.fig = add_hinge(self.fig, HEIGHT, self._ev(link['hinge']))
            if 'mobile_support' in link:
                self.fig = add_mobile_support(self.fig, HEIGHT, self._ev(link['mobile_support']))
            if 'fixed_support' in link:
                self.fig = add_fixed_support(self.fig, HEIGHT, self._ev(link['fixed_support']))
            if 'roller' in link:
                self.fig = add_roller_support(self.fig, HEIGHT, self._ev(link['roller']))                
        
        # Add shear forces
        for force in self.shear_forces:
            if force['n'] < 0:
                if self._ev(force['start']) < self.beam_size:
                    if force['pos']:
                        if force['r']:
                            self.fig = add_vector(self.fig, (self._ev(force['start']), -1.5*HEIGHT + 0.15*HEIGHT), (self._ev(force['start']), -HEIGHT/2 + 0.15*HEIGHT), format_subs(force['value']), 'red')
                        else:
                            self.fig = add_vector(self.fig, (self._ev(force['start']), HEIGHT), (self._ev(force['start']), 2*HEIGHT), format_subs(force['value']))
                    else:
                        if force['r']:
                            self.fig = add_vector(self.fig, (self._ev(force['start']), HEIGHT/2), (self._ev(force['start']), -1.5*HEIGHT), format_subs(force['value']), 'red')
                        else:
                            self.fig = add_vector(self.fig, (self._ev(force['start']), 2.45*HEIGHT), (self._ev(force['start']), HEIGHT), format_subs(force['value']))
                else:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), 2.45*HEIGHT), (self._ev(force['start']), HEIGHT), format_subs(force['value']))                    
                    else:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), HEIGHT), (self._ev(force['start']), 2*HEIGHT), format_subs(force['value']))
            else:
                x = np.linspace(self._ev(force['start']), self._ev(force['stop']), 100)
                if force['n'] == 0:
                    y = [2* HEIGHT for i in range(len(x))]
                else:
                    y = (x - self._ev(force['start'])) ** force['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, HEIGHT, format_subs(force['value']), up=force['pos'])
        
        # Add normal forces
        for force in self.normal_forces:
            if force['n'] < 0:
                if self._ev(force['start']) > 0:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), HEIGHT/2), (self._ev(force['start']) + HEIGHT, HEIGHT/2), format_subs(force['value']))
                    else:
                        self.fig = add_vector(self.fig, (self._ev(force['start']) + HEIGHT, HEIGHT/2), (self._ev(force['start']), HEIGHT/2), format_subs(force['value']))
                else:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (self._ev(force['start']) + HEIGHT, HEIGHT/2), (self._ev(force['start']), HEIGHT/2), format_subs(force['value']))                 
                    else:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), HEIGHT/2), (self._ev(force['start']) + HEIGHT, HEIGHT/2), format_subs(force['value']))
            else:
                x = np.linspace(self._ev(force['start']), self._ev(force['stop']), 100)
                if force['n'] == 0:
                    y = [2* HEIGHT for i in range(len(x))]
                else:
                    y = (x - self._ev(force['start'])) ** force['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, HEIGHT, format_subs(force['value']), up=force['pos'], side=True)
        
        # Add twisting moments
        for moment in self.twisting_moments:
            if moment['n'] < 0:
                if self._ev(moment['start']) > 0:
                    if moment['pos']:
                        self.fig = add_vector(self.fig, (self._ev(moment['start']), HEIGHT/2), (self._ev(moment['start']) + 0.75*HEIGHT, HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 0.7*HEIGHT, HEIGHT/2), (self._ev(moment['start']) + 1*HEIGHT, HEIGHT/2))
                    else:
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 0.75*HEIGHT, HEIGHT/2), (self._ev(moment['start']), HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 1*HEIGHT, HEIGHT/2), (self._ev(moment['start']) + 0.7*HEIGHT, HEIGHT/2))
                else:
                    if moment['pos']:
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 0.75*HEIGHT, HEIGHT/2), (self._ev(moment['start']), HEIGHT/2), format_subs(moment['value'])) 
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 1*HEIGHT, HEIGHT/2), (self._ev(moment['start']) + 0.7*HEIGHT, HEIGHT/2))                
                    else:
                        self.fig = add_vector(self.fig, (self._ev(moment['start']), HEIGHT/2), (self._ev(moment['start']) + 0.75*HEIGHT, HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 0.7*HEIGHT, HEIGHT/2), (self._ev(moment['start']) + 1*HEIGHT, HEIGHT/2))
            else:
                x = np.linspace(self._ev(moment['start']), self._ev(moment['stop']), 100)
                if moment['n'] == 0:
                    y = [2* HEIGHT for i in range(len(x))]
                else:
                    y = (x - self._ev(moment['start'])) ** moment['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, HEIGHT, format_subs(moment['value']), up=moment['pos'], side=True, double=True)    

        # Add bending moments
        for moment in self.bending_moments:
            if self._ev(moment['start']) < self.beam_size:
                if moment['pos']:
                    self.fig = add_semicircle_arrow(self.fig, self._ev(moment['start']), HEIGHT/2, orientation='clockwise', label=format_subs(moment['value']))
                else:
                    self.fig = add_semicircle_arrow(self.fig, self._ev(moment['start']), HEIGHT/2, label=format_subs(moment['value']))
            else:
                if moment['pos']:
                    self.fig = add_semicircle_arrow(self.fig, self._ev(moment['start']), HEIGHT/2, label=format_subs(moment['value']))                 
                else:
                    self.fig = add_semicircle_arrow(self.fig, self._ev(moment['start']), HEIGHT/2, orientation='clockwise', label=format_subs(moment['value']))

        return self.fig