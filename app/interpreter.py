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
        "L": 5
    },
    "points": {
        "A": 0,
        "B": 'L'
    },
    "links": [],
    "shear_forces": [],
    "normal_forces": [],
    "twisting_moments": [],
    "bending_moments": [],
    "boundary_conditions": []
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
        self.beam_size: float = b['beam_size']
        self.variables: dict = b['variables']
        self.points: dict = b['points']
        self.links: list = b['links']
        self.shear_forces: list = b['shear_forces']
        self.normal_forces: list = b['normal_forces']
        self.twisting_moments: list = b['twisting_moments']
        self.bending_moments: list = b['bending_moments']
        self.boundary_conditions: list = b['boundary_conditions']
        self.protected_symbols = [
            'N_x',
            'M_x',
            'M_z',
            'V_y',
            'u',
            'q',
            'p',
            't',
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
        b['boundary_conditions'] = self.boundary_conditions
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
                't',
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
                self.add_shear_force(f'R_y_{link_point}', 0, position, position, -1, True, True)
                self._enforce_reaction_consistency()
            elif (0 < self._ev(position) < self.beam_size):
                link_point = self._get_and_add_point(position)
                self._enforce_reaction_consistency()
            self.links.append({link_type: position})
            return True
        return False

    def _create_load(self, value: str, value_min: str, start: str, stop: str, n: int, pos: bool = True, reaction: bool = False) -> dict:
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
        if value_min == '':
            value_min = 0
        if value_min not in self.variables:
            try:
                float(value_min)
            except:
                if reaction:
                    if reaction not in self.protected_symbols:
                        self.protected_symbols.append(value_min)
                else:
                    value_min = self._update_variables(value_min)                    
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
                "value_min": value_min,
                "start": start,
                "stop": stop,
                "n": n,
                "pos": pos,
                "r": reaction
            }          

    def add_shear_force(self, value: str, value_min: str, start: str, stop: str, n: int, pos: bool = True, reaction: bool = False):
        """
        Adds a concentrated or distributed shear force.   

        Parameters:
            value        : Value of the force.        
            start         : Start position of the domain of the force.
            stop         : End position of the domain of the force.
            n         : Exponent of the polynomial representing the force.
            pos         : True if the force is positive, False otherwise.
        """ 
        if 'I_zz' not in self.variables or 'E' not in self.variables:
            self.add_variable('I_zz',1)
            self.add_variable('E',200 * 1e9)        
        self.shear_forces.append(self._create_load(value, value_min, start, stop, n, pos, reaction))  

    def add_normal_force(self, value: str, value_min: str, start: str, stop: str, n: int, pos: bool = True, reaction: bool = False):
        """
        Adds a concentrated or distributed normal force.   

        Parameters:
            value        : Value of the force.        
            start         : Start position of the domain of the force.
            stop         : End position of the domain of the force.
            n         : Exponent of the polynomial representing the force.
            pos         : True if the force is positive, False otherwise.
        """       
        if 'A' not in self.variables or 'E' not in self.variables:
            self.add_variable('A',1)
            self.add_variable('E',200 * 1e9)        
        self.normal_forces.append(self._create_load(value, value_min, start, stop, n, pos, reaction))  

    def add_twisting_moment(self, value: str, value_min: str, start: str, stop: str, n: int, pos: bool = True, reaction: bool = False):
        """
        Adds a concentrated or distributed twisting moment.   

        Parameters:
            value        : Value of the moment.        
            start         : Start position of the domain of the moment.
            stop         : End position of the domain of the moment.
            n         : Exponent of the polynomial representing the moment.
            pos         : True if the moment is positive, False otherwise.
        """        
        if 'J_p' not in self.variables or 'G' not in self.variables:
            self.add_variable('J_p',1)
            self.add_variable('G', 76 * 1e9)        
        self.twisting_moments.append(self._create_load(value, value_min, start, stop, n, pos, reaction))  

    def add_bending_moment(self, value: str, value_min: str, start: str, stop: str, n: int, pos: bool = True, reaction: bool = False):
        """
        Adds a concentrated or distributed bending moment.   

        Parameters:
            value        : Value of the moment.        
            start         : Start position of the domain of the moment.
            stop         : End position of the domain of the moment.
            n         : Exponent of the polynomial representing the moment.
            pos         : True if the moment is positive, False otherwise.
        """   
        if 'I_zz' not in self.variables or 'E' not in self.variables:
            self.add_variable('I_zz',1)
            self.add_variable('E',200 * 1e9)                
        self.bending_moments.append(self._create_load(value, value_min, start, stop, n, pos, reaction))
    
    def _get_load_equation(self, force: dict):
        x = symbols('x')
        try:
            p = float(force['value'])
        except:
            p = self._s(force['value'])
        try:
            q = float(force['value_min'])
        except:
            q = self._s(force['value_min'])        
        if not force['pos']:
            p = -p
            q = -q
        x_ = self._s(f"({force['stop']})-({force['start']})")
        y_ = (p-q)
        r = y_/x_
        if self._ev(force['stop']) < self.beam_size and force['n'] >= 0:
            if force['n'] <= 0:
                eq = p * SingularityFunction(x, self._s(force['start']), force['n']) - p * SingularityFunction(x, self._s(force['stop']), force['n'])
            else:
                if self._ev(force['start']) == 0 and q == 0:
                    eq = r * SingularityFunction(x, self._s(force['start']), force['n']) - r * SingularityFunction(x, self._s(force['stop']), force['n'])
                else:
                    eq = q * SingularityFunction(x, self._s(force['start']), 0) + r * SingularityFunction(x, self._s(force['start']), force['n']) - (q * SingularityFunction(x, self._s(force['stop']), 0) + r * SingularityFunction(x, self._s(force['stop']), force['n']))
        else:
            if force['n'] <= 0:
                eq = p * SingularityFunction(x, self._s(force['start']), force['n'])
            else:
                if self._ev(force['start']) == 0 and q == 0:
                    eq = r * SingularityFunction(x, self._s(force['start']), force['n'])
                else:
                    eq = q * SingularityFunction(x, self._s(force['start']), 0) + r * SingularityFunction(x, self._s(force['start']), force['n'])
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
            if load['n'] < 0 and (self._ev(load['start']) == 0 or self._ev(load['start']) == self.beam_size):
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
        s, f = self._get_functions_and_symbols()
        return [latex_with_threshold(Eq(f[eq[0]](self._s(eq[1])),self._s(eq[2]))) for eq in self.boundary_conditions]

    def _get_position_constants(self) -> List[str]:
        """
        """
        consts = set(['L'])
        for variable in self.variables:
            for load in self.normal_forces + self.twisting_moments + self.shear_forces + self.bending_moments:
                if isinstance(load['start'],str):
                    if variable in load['start']:
                        consts.add(variable)
                if isinstance(load['stop'],str):
                    if variable in load['stop']:
                        consts.add(variable)
            for link in self.links:
                if isinstance(link[list(link.keys())[0]], str):
                    if variable in link[list(link.keys())[0]]:
                        consts.add(variable)
        return list(consts)

    def _check_boundaries(self, left=True, kind='bending') -> bool:
        """
        """
        lds = {
            'normal': self.normal_forces,
            'twisting': self.twisting_moments,
            'shear': self.shear_forces,
            'bending': self.bending_moments
        }
        t = 0 if left else self.beam_size
        for load in lds[kind]:
            if self._ev(load['start']) == t and load['n'] < 0:
                return True
        for link in self.links:
            if self._ev(link[list(link.keys())[0]]) == t:
                return True            
        return False

    def calculate_boundary_conditions(self):
        """
        """ 
        remove_zero_flags(self.boundary_conditions)
        s, f = self._get_functions_and_symbols()
        mapped_points = {name: [] for name in [str(fun) for fun in f]}

        for item in self.boundary_conditions:
            mapped_points[item[0]].append(item[1])

        degrees_of_freedom = {
            'shear': 4,
            'normal': 2,
            'twisting': 2
        }

        def add_condition(func: str, pos: str, value: str) -> None:
            """
            Adds a boundary condition   

            Parameters:
                func        : String of the function symbol.
                pos        : String of the position.

            Returns:
                out: None          
            """
            add_item(self.boundary_conditions, [func, pos, value, 0])
            mapped_points[func].append(pos) if pos not in mapped_points[func] else None

        def create_boundary_condition_from_force(force: dict, func: str, kind: str) -> None:
            if force['r']:
                degrees_of_freedom[kind] += 1
            if force['n'] < 0:
                if self._ev(force['start']) == 0:
                    if force['pos'] and 0 not in mapped_points[func]:
                        add_condition(func, 0, force['value'])
                    elif 0 not in mapped_points[func]:
                        add_condition(func, 0, f"-({force['value']})")
                elif self._ev(force['start']) == self.beam_size:
                    if force['pos'] and 'L' not in mapped_points[func]:
                        add_condition(func, 'L', force['value'])
                    elif 'L' not in mapped_points[func]:
                        add_condition(func, 'L', f"-({force['value']})")

        # Normal
        for force in self.normal_forces:
            create_boundary_condition_from_force(force, 'N_x', 'normal')

        # Twisting
        for moment in self.twisting_moments:
            create_boundary_condition_from_force(moment, 'M_x', 'twisting')

        # Shear
        for force in self.shear_forces:
            create_boundary_condition_from_force(force, 'V_y', 'shear')

        # Bending
        for moment in self.bending_moments:
            create_boundary_condition_from_force(moment, 'M_z', 'shear')

        # Links
        for link in self.links:
            link_type = list(link.keys())[0]
            position = link[link_type]
            if link_type == 'cantilever':
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    # Shear and bending
                    add_condition('theta_Z', position, 0)
                    add_condition('v', position, 0)    
                if len(self.normal_forces) > 0:
                    # Normal     
                    add_condition('u', position, 0)
                if len(self.twisting_moments) > 0:
                    # Twisting                           
                    add_condition('phi', position, 0)                        
            if link_type == 'mobile_support':
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    # Shear and bending
                    if self._ev(position) == 0 or self._ev(position) == self.beam_size:
                        add_condition('M_z', position, 0)
                    add_condition('v', position, 0)
            if link_type == 'fixed_support':
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    # Shear and bending
                    if self._ev(position) == 0 or self._ev(position) == self.beam_size:
                        add_condition('M_z', position, 0)
                    add_condition('v', position, 0)
                if len(self.normal_forces) > 0:
                    # Normal 
                    add_condition('u', position, 0)           
            if link_type == 'hinge':
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    add_condition('M_z', position, 0)   
            if link_type == 'roller':
                if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                    add_condition('V_y', position, 0)

        # Empty boundaries
        if len(mapped_points['M_z']) + len(mapped_points['V_y']) + len(mapped_points['theta_Z']) + len(mapped_points['v']) < degrees_of_freedom['shear']:
            if len(self.shear_forces) > 0 or len(self.bending_moments) > 0:
                if 0 not in mapped_points['V_y'] and not self._check_boundaries(True,'shear'):
                    add_condition('V_y', 0, 0)
                if 'L' not in mapped_points['V_y'] and not self._check_boundaries(False,'shear'):
                    add_condition('V_y', 'L', 0)           
                if 0 not in mapped_points['M_z'] and not self._check_boundaries():
                    add_condition('M_z', 0, 0)
                if 'L' not in mapped_points['M_z'] and not self._check_boundaries(False):
                    add_condition('M_z', 'L', 0)       

        if len(mapped_points['N_x']) + len(mapped_points['u']) < degrees_of_freedom['normal']:
            if len(self.normal_forces) > 0:
                if 0 not in mapped_points['N_x'] and not self._check_boundaries(True,'normal'):
                    add_condition('N_x', 0, 0)
                if 'L' not in mapped_points['N_x'] and not self._check_boundaries(False,'normal'):
                    add_condition('N_x', 'L', 0)
                       
        if len(mapped_points['M_x']) + len(mapped_points['phi']) < degrees_of_freedom['twisting']:
            if len(self.twisting_moments) > 0:
                if 0 not in mapped_points['M_x'] and not self._check_boundaries(True,'twisting'):
                    add_condition('M_x', 0, 0)
                if 'L' not in mapped_points['M_x'] and not self._check_boundaries(False,'twisting'):
                    add_condition('M_x', 'L', 0)                       

        priority = ['v', 'theta_Z', 'M_z', 'V_y', 'u', 'N_x', 'phi', 'M_x']

        # Build a quick lookup dict so that each name maps to its index (0..7)
        priority_index = {name: i for i, name in enumerate(priority)}

        # Sort so that:
        #    - Items with item[1] == 0 come first (because (item[1] != 0) is False → 0)
        #    - Within the zero‐group and the nonzero‐group, they follow the order in `priority`
        self.boundary_conditions.sort(
            key=lambda item: (
                # first key: False (0) if item[1]==0, True (1) otherwise
                item[1] != 0,
                # second key: the index in our priority list, or a large default if not found
                priority_index.get(item[0], len(priority))
            )
        )

        enough = {
            'shear': len(mapped_points['M_z']) + len(mapped_points['V_y']) + len(mapped_points['theta_Z']) + len(mapped_points['v']) >= degrees_of_freedom['shear'],
            'normal': len(mapped_points['N_x']) + len(mapped_points['u']) >= degrees_of_freedom['normal'],
            'twisting': len(mapped_points['M_x']) + len(mapped_points['phi']) >= degrees_of_freedom['twisting']
        }

        return mapped_points, enough

    def solve(self):
        """
        """
        # WARNING: insanely ugly code ahead
        s, f = self._get_functions_and_symbols()
        mp, en = self.calculate_boundary_conditions()
        x = symbols('x')
        vardict = self._get_symbol_value()
        relations = []
        # force position symbols to be positive
        positive_symbols = set()
        for load in self.normal_forces + self.shear_forces + self.bending_moments + self.twisting_moments:
            positive_symbols |= self._s(load['start']).free_symbols | self._s(load['stop']).free_symbols
        for link in self.links:
            positive_symbols |= self._s(link[list(link.keys())[0]]).free_symbols
        for symb in positive_symbols:
            sb = Symbol(str(symb), real=True, positive=True)
            vardict[sb] = vardict.pop(self.symbols[str(symb)])
        for item in vardict:
            if str(item) in self.variables:
                for other in vardict:
                    if str(other) in self.variables:
                        if self.variables[str(item)] > self.variables[str(other)]:
                            relations.append(Gt(item, other))
        symdict = {str(sb): sb for sb in vardict}
        solution_blocks = []

        posdict = {}
        posconst = self._get_position_constants()
        for item in posconst:
            posdict[symdict[item]] = vardict[symdict[item]]

        print(posdict)

        print(self.boundary_conditions, '\n')
        print(mp, en, '\n')

        def calculate_constant(cond: Item, constant_det: dict, constant_value: dict, p: Expr, k: Expr, constant_lit: dict):
            # sb = position of condition
            if isinstance(cond[1], str):
                sb = sympify(cond[1], locals=symdict)
            else:
                sb = cond[1]
            # isolate the constant
            sol = solve(Eq(sympify(cond[2], locals=symdict), p),k)[0]
            constant_det[str(k)].append(latex_with_threshold(collapse_singularity_functions(Eq(k, sol).subs(x, sb), x, relations)))  
            expr = Eq(k, sol).subs(x, sb)
            # substitute positions
            for v in posdict:
                expr = expr.subs(v, posdict[v])            
            constant_det[str(k)].append(latex_with_threshold(collapse_singularity_functions(expr, x, relations)))
            constant_lit[str(k)] = simplify(expr.rhs)
            # substitute known constants
            for const in constant_value:
                sol = sol.subs(const, constant_value[const])            
            # evaluate
            if isinstance(sb, Expr):
                sol = Eq(k, sol).subs(x, sb)
                for s in sb.free_symbols:
                    sol = sol.subs(s, vardict[s]).evalf()
            else:
                sol = Eq(k, sol).subs(x, sb).evalf()
            # replace symbols not in the position string with their numerical value
            for v in vardict:
                if str(v) != str(sb):
                    sol = sol.subs(v, vardict[v], evaluate=False)
            # isolate again if needed
            if k in sol.rhs.free_symbols:
                sol = Eq(k, solve(sol, k)[0])
            constant_det[str(k)].append(latex_with_threshold(sol))
            constant_value[str(k)] = sol.rhs.evalf(subs=vardict)
            # remove duplicates
            constant_det[str(k)] = list(dict.fromkeys(constant_det[str(k)])) 

        def finalize_constants(constant_strings: list, reaction_strings: list, reactions: list, constant_value: dict, constant_det: dict):
            for k in constant_strings:
                if k in constant_value:
                    try: 
                        float(constant_value[str(k)])
                    except:
                        sol = constant_value[str(k)]
                        for const in constant_value:
                            if const != k:
                                sol = sol.subs(const, constant_value[const])
                        constant_value[str(k)] = sol.evalf(subs=vardict)
                        final_constant_expression = latex_with_threshold(Eq(c[constant_strings.index(str(k))],constant_value[k]))
                        if final_constant_expression not in constant_det[k]:
                            constant_det[k].append(final_constant_expression)
            for k in reaction_strings:
                if k in constant_value:
                    try: 
                        float(constant_value[str(k)])
                    except:             
                        sol = constant_value[str(k)]
                        for const in constant_value:
                            if const != k:
                                sol = sol.subs(const, constant_value[const])
                        constant_value[str(k)] = sol.evalf(subs=vardict)                                                       
                        final_constant_expression = latex_with_threshold(Eq(reactions[reaction_strings.index(str(k))],constant_value[k]))
                        if final_constant_expression not in constant_det[k]:
                            constant_det[k].append(final_constant_expression)            

        def plot_equations(plots: dict, final_eqs: list, constants: list, constant_value: dict, constant_lit: dict):
            for i in range(len(final_eqs)):
                x_axis = np.linspace(0, self.variables['L'],1000)
                eq = final_eqs[i]
                eqlit = copy(final_eqs[i])
                for k in constants:
                    eq = eq.subs(k,constant_value[str(k)])
                    eqlit = eqlit.subs(k, constant_lit[str(k)])
                for r in reactions:
                    eq = eq.subs(r,constant_value[str(r)])
                    eqlit = eqlit.subs(k, constant_lit[str(r)])
                eqlit = collapse_singularity_functions(eqlit, x, relations)
                eq = eq.evalf(subs=self._remove_protected_symbols(vardict))
                f = lambdify(x, eq.rhs, modules=[mapping, 'numpy'])
                y_axis = f(x_axis)
                if isinstance(y_axis, float):
                    y_axis = np.zeros(len(x_axis)) + y_axis
                title = format_label(str(eq.lhs)).replace('*','')
                plots[title] = fig_to_rotated_img(create_filled_line_figure(
                    x_axis,
                    y_axis,
                    title
                ))
                final_eqs[i] = latex_with_threshold(eqlit)
                final_eqs.append(latex_with_threshold(eq))

        # Normal
        if len(self.normal_forces) > 0 and en['normal']:
            normal_block = {'name': 'Força Normal'}
            reactions = []
            # create load equation
            p = sympify(0)
            for force in self.normal_forces:
                if (force['n'] < 0 and (0 < self._ev(force['start']) < self.beam_size)) or force['n'] >= 0:
                    p += self._get_load_equation(force)
                if force['r']:
                    reactions.append(unify_symbols(self._s(force['value']), symdict))
            p = unify_symbols(p, symdict)
            normal_block['load_equation'] = latex_with_threshold(Eq(f['p'](x), p))

            c = []
            ps = []            
            if len(mp['u']) > 0:
                normal_block['differential_equation'] = latex_with_threshold(Eq(Eq(s['A']*s['E'] * Derivative(f['u'](x),(x,2)), -f['p'](x)), -p, evaluate=False))
                steps = []
                final_eqs = []
                plots = {}                
                c.append(Symbol('c_1')) # first constant
                p1 = integrate(-p, x) + c[0] # first integration
                ps.append(p1)
                steps.append(latex_with_threshold(Eq(Eq(s['A']*s['E'] * Derivative(f['u'](x),(x,1)),f['N_x'](x)), p1,evaluate=False)))
                final_eqs.append(Eq(f['N_x'](x), p1, evaluate=False))
                c.append(Symbol('c_2')) # second constant 
                p2 = integrate(p1, x) + c[1] # second integration
                ps.append(p2)
                steps.append(latex_with_threshold(Eq(s['A']*s['E'] * f['u'](x),p2,evaluate=False)))
                final_eqs.append(Eq(s['A']*s['E'] * f['u'](x),p2,evaluate=False))
                normal_block['integration_steps'] = steps
                # find out constants
                constant_det = {}
                constant_value = {}
                constant_lit = {}
                total_unknowns = len(c) + len(reactions)
                constant_strings = [str(g) for g in c]
                reaction_strings = [str(g) for g in reactions]
                for cond in self.boundary_conditions:                    
                    if cond[0] == 'N_x':
                        k = c[0]
                        p = ps[0]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[:1]):
                                if item not in constant_det and item != ogk:                                
                                    k = c[i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings[:1]) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                                  
                        constant_det[str(k)] = []
                        # Nx(position) = value = first integration of -p(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                    if cond[0] == 'u':
                        k = c[1]
                        p = (1/(s['A']*s['E'])) * ps[1]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings):
                                if item not in constant_det and item != ogk:                                
                                    k = c[i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                       
                        constant_det[str(k)] = []
                        # u(position) = value = second integration of -p(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break       
                # finalize constants
                finalize_constants(constant_strings, reaction_strings, reactions, constant_value, constant_det)
                normal_block['constants'] = constant_det
                # plots
                plot_equations(plots, final_eqs, c, constant_value, constant_lit)
                normal_block['final_equations'] = final_eqs
                normal_block['plots'] = plots
            elif len(mp['N_x']) > 0:
                normal_block['differential_equation'] = latex_with_threshold(Eq(Eq(Derivative(f['N_x'](x),(x,1)), -f['p'](x)), -p, evaluate=False))
                steps = []
                final_eqs = []
                plots = {}                
                c.append(Symbol('c_1')) # first constant
                p1 = integrate(-p, x) + c[0] # first integration
                ps.append(p1)
                steps.append(latex_with_threshold(Eq(f['N_x'](x), p1,evaluate=False)))
                final_eqs.append(Eq(f['N_x'](x), p1, evaluate=False))
                normal_block['integration_steps'] = steps
                # find out constants
                constant_det = {}
                constant_value = {}
                constant_lit = {}
                total_unknowns = len(c) + len(reactions)
                constant_strings = [str(g) for g in c]
                reaction_strings = [str(g) for g in reactions]
                for cond in self.boundary_conditions:
                    if cond[0] == 'N_x':
                        k = c[0]
                        p = ps[0]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[:1]):
                                if item not in constant_det and item != ogk:                                
                                    k = c[i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings[:1]) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                                  
                        constant_det[str(k)] = []
                        # Nx(position) = value = first integration of -p(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                # finalize constants
                finalize_constants(constant_strings, reaction_strings, reactions, constant_value, constant_det)
                normal_block['constants'] = constant_det
                # plots
                plot_equations(plots, final_eqs, c, constant_value, constant_lit)
                normal_block['final_equations'] = final_eqs
                normal_block['plots'] = plots                                           

            solution_blocks.append(normal_block)

        # Twisting
        if len(self.twisting_moments) > 0 and en['twisting']:
            twisting_block = {'name': 'Torção de Seções Circulares'}
            reactions = []
            # create load equation
            p = sympify(0)
            for force in self.twisting_moments:
                if (force['n'] < 0 and (0 < self._ev(force['start']) < self.beam_size)) or force['n'] >= 0:
                    p += self._get_load_equation(force)
                if force['r']:
                    reactions.append(unify_symbols(self._s(force['value']), symdict))
            p = unify_symbols(p, symdict)
            twisting_block['load_equation'] = latex_with_threshold(Eq(f['t'](x), p))

            c = []
            ps = []            
            if len(mp['phi']) > 0:
                twisting_block['differential_equation'] = latex_with_threshold(Eq(Eq(s['J_p']*s['G'] * Derivative(f['phi'](x),(x,2)), -f['t'](x)), -p, evaluate=False))
                steps = []
                final_eqs = []
                plots = {}                
                c.append(Symbol('c_1')) # first constant
                p1 = integrate(-p, x) + c[0] # first integration
                ps.append(p1)
                steps.append(latex_with_threshold(Eq(Eq(s['J_p']*s['G'] * Derivative(f['phi'](x),(x,1)),f['M_x'](x)), p1,evaluate=False)))
                final_eqs.append(Eq(f['M_x'](x), p1, evaluate=False))
                c.append(Symbol('c_2')) # second constant 
                p2 = integrate(p1, x) + c[1] # second integration
                ps.append(p2)
                steps.append(latex_with_threshold(Eq(s['J_p']*s['G'] * f['phi'](x),p2,evaluate=False)))
                final_eqs.append(Eq(s['J_p']*s['G'] * f['phi'](x),p2,evaluate=False))
                twisting_block['integration_steps'] = steps
                # find out constants
                constant_det = {}
                constant_value = {}
                constant_lit = {}
                total_unknowns = len(c) + len(reactions)
                constant_strings = [str(g) for g in c]
                reaction_strings = [str(g) for g in reactions]
                for cond in self.boundary_conditions:
                    if cond[0] == 'M_x':
                        k = c[0]
                        p = ps[0]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[:1]):
                                if item not in constant_det and item != ogk:                                
                                    k = c[i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings[:1]) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                                  
                        constant_det[str(k)] = []
                        # Mx(position) = value = first integration of -t(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                    if cond[0] == 'phi':
                        k = c[1]
                        p = (1/(s['J_p']*s['G'])) * ps[1]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings):
                                if item not in constant_det and item != ogk:                                
                                    k = c[i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                         
                        constant_det[str(k)] = []
                        # phi(position) = value = second integration of -t(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break       
                # finalize constants
                finalize_constants(constant_strings, reaction_strings, reactions, constant_value, constant_det)
                twisting_block['constants'] = constant_det
                # plots
                plot_equations(plots, final_eqs, c, constant_value, constant_lit)
                twisting_block['final_equations'] = final_eqs
                twisting_block['plots'] = plots
            elif len(mp['M_x']) > 0:
                twisting_block['differential_equation'] = latex_with_threshold(Eq(Eq(Derivative(f['M_x'](x),(x,1)), -f['t'](x)), -p, evaluate=False))
                steps = []
                final_eqs = []
                plots = {}                
                c.append(Symbol('c_1')) # first constant
                p1 = integrate(-p, x) + c[0] # first integration
                ps.append(p1)
                steps.append(latex_with_threshold(Eq(f['M_x'](x), p1,evaluate=False)))
                final_eqs.append(Eq(f['M_x'](x), p1, evaluate=False))
                twisting_block['integration_steps'] = steps
                # find out constants
                constant_det = {}
                constant_value = {}
                constant_lit = {}
                total_unknowns = len(c) + len(reactions)
                constant_strings = [str(g) for g in c]
                reaction_strings = [str(g) for g in reactions]
                for cond in self.boundary_conditions:
                    if cond[0] == 'M_x':
                        k = c[0]
                        p = ps[0]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[:1]):
                                if item not in constant_det and item != ogk:                                
                                    k = c[i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings[:1]) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                                 
                        constant_det[str(k)] = []
                        # Mx(position) = value = first integration of -t(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                # finalize constants
                finalize_constants(constant_strings, reaction_strings, reactions, constant_value, constant_det)
                twisting_block['constants'] = constant_det
                # plots
                plot_equations(plots, final_eqs, c, constant_value, constant_lit)
                twisting_block['final_equations'] = final_eqs
                twisting_block['plots'] = plots                                           

            solution_blocks.append(twisting_block)            

        # Shear and bending
        if len(self.shear_forces) > 0 or len(self.bending_moments) > 0 and en['shear']:
            shear_block = {'name': 'Flexão Pura'}
            reactions = []
            # create load equation
            p = sympify(0)
            for force in self.shear_forces + self.bending_moments:
                if (force['n'] < 0 and (0 < self._ev(force['start']) < self.beam_size)) or force['n'] >= 0:
                    p += self._get_load_equation(force)
                if force['r']:
                    reactions.append(unify_symbols(self._s(force['value']), symdict))
            p = unify_symbols(p, symdict)
            shear_block['load_equation'] = latex_with_threshold(Eq(f['q'](x), p))

            c = []
            ps = []            
            if len(mp['v']) > 0:
                shear_block['differential_equation'] = latex_with_threshold(Eq(Eq(s['I_zz']*s['E'] * Derivative(f['v'](x),(x,4)), f['q'](x)), p, evaluate=False))
                steps = []
                final_eqs = []
                plots = {}                
                c.append(Symbol('c_1')) # first constant
                p1 = integrate(p, x) + c[0] # first integration
                ps.append(p1)
                steps.append(latex_with_threshold(Eq(Eq(s['I_zz']*s['E'] * Derivative(f['v'](x),(x,3)),f['V_y'](x)), p1,evaluate=False)))
                final_eqs.append(Eq(f['V_y'](x), p1, evaluate=False))
                c.append(Symbol('c_2')) # second constant 
                p2 = integrate(p1, x) + c[1] # second integration
                ps.append(p2)
                steps.append(latex_with_threshold(Eq(Eq(s['I_zz']*s['E'] * Derivative(f['v'](x),(x,2)),f['M_z'](x)), p2,evaluate=False)))
                final_eqs.append(Eq(f['M_z'](x), p2, evaluate=False))
                c.append(Symbol('c_3')) # third constant 
                p3 = integrate(p2, x) + c[2] # third integration
                ps.append(p3)
                steps.append(latex_with_threshold(Eq(Eq(s['I_zz']*s['E'] * Derivative(f['v'](x),(x,1)),s['I_zz']*s['E'] * f['theta_Z'](x)), p3,evaluate=False)))
                final_eqs.append(Eq(s['I_zz']*s['E'] * f['theta_Z'](x), p3, evaluate=False))
                c.append(Symbol('c_4')) # fourth constant 
                p4 = integrate(p3, x) + c[3] # fourth integration
                ps.append(p4)
                steps.append(latex_with_threshold(Eq(s['I_zz']*s['E'] * f['v'](x), p4,evaluate=False)))
                final_eqs.append(Eq(s['I_zz']*s['E'] * f['v'](x), p4, evaluate=False))                                   
                shear_block['integration_steps'] = steps
                # find out constants
                constant_det = {}
                constant_value = {}
                constant_lit = {}
                total_unknowns = len(c) + len(reactions)
                constant_strings = [str(g) for g in c]
                reaction_strings = [str(g) for g in reactions]
                for cond in self.boundary_conditions:
                    if cond[0] == 'V_y':
                        k = c[0]
                        p = ps[0]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[:1][::-1]):
                                if item not in constant_det and item != ogk:                                
                                    k = c[:1][::-1][i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings[:1]) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                                  
                        constant_det[str(k)] = []
                        # V_y(position) = value = first integration of q(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                    if cond[0] == 'M_z':
                        k = c[1]
                        p = ps[1]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[:2][::-1]):                                
                                if item not in constant_det and item != ogk:                                
                                    k = c[:2][::-1][i]
                                    break   
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings[:2]) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                     
                        constant_det[str(k)] = []
                        # M_z(position) = value = second integration of q(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break  
                    if cond[0] == 'theta_Z':
                        k = c[2]
                        p = ps[2]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[:3][::-1]):
                                if item not in constant_det and item != ogk:                                
                                    k = c[:3][::-1][i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings[:3]) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                  
                        constant_det[str(k)] = []
                        # theta_Z(position) = value = third integration of q(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                    if cond[0] == 'v':
                        k = c[3]
                        p = ps[3]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[::-1]):                                
                                if item not in constant_det and item != ogk:                                
                                    k = c[::-1][i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                         
                        constant_det[str(k)] = []
                        # v(position) = value = fourth integration of q(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break                                                            
                # finalize constants
                finalize_constants(constant_strings, reaction_strings, reactions, constant_value, constant_det)
                shear_block['constants'] = constant_det
                # plots
                plot_equations(plots, final_eqs, c, constant_value, constant_lit)
                shear_block['final_equations'] = final_eqs
                shear_block['plots'] = plots
            elif len(mp['M_z']) > 0:
                shear_block['differential_equation'] = latex_with_threshold(Eq(Eq(Derivative(f['M_z'](x),(x,2)), f['q'](x)), p, evaluate=False))
                steps = []
                final_eqs = []
                plots = {}                
                c.append(Symbol('c_1')) # first constant
                p1 = integrate(p, x) + c[0] # first integration
                ps.append(p1)
                steps.append(latex_with_threshold(Eq(Eq(Derivative(f['v'](x),(x,1)),f['V_y'](x)), p1,evaluate=False)))
                final_eqs.append(Eq(f['V_y'](x), p1, evaluate=False))
                c.append(Symbol('c_2')) # second constant 
                p2 = integrate(p1, x) + c[1] # second integration
                ps.append(p2)
                steps.append(latex_with_threshold(Eq(f['M_z'](x), p2,evaluate=False)))
                final_eqs.append(Eq(f['M_z'](x), p2, evaluate=False))                                
                shear_block['integration_steps'] = steps
                # find out constants
                constant_det = {}
                constant_value = {}
                constant_lit = {}
                total_unknowns = len(c) + len(reactions)
                constant_strings = [str(g) for g in c]
                reaction_strings = [str(g) for g in reactions]
                for cond in self.boundary_conditions:
                    if cond[0] == 'V_y':
                        k = c[0]
                        p = ps[0]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings[:1]):
                                if item not in constant_det and item != ogk:                                
                                    k = c[i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings[:1]) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]
                                    break                                
                        constant_det[str(k)] = []
                        # V_y(position) = value = first integration of q(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break
                    if cond[0] == 'M_z':
                        k = c[1]
                        p = ps[1]
                        ogk = str(k)
                        # check if constant has already been determined, if true then get next one
                        if ogk in constant_det:
                            for i, item in enumerate(constant_strings):
                                if item not in constant_det and item != ogk:                                
                                    k = c[i]
                                    break
                        # if constants have been determined, find out reactions
                        if str(k) in constant_det and constant_strings.index(str(k)) == len(constant_strings) - 1:
                            for i, item in enumerate(reaction_strings):
                                if item not in constant_det:
                                    k = reactions[i]                                    
                                    break                      
                        constant_det[str(k)] = []
                        # M_z(position) = value = second integration of q(x)
                        constant_det[str(k)].append(latex_with_threshold(Eq(Eq(f[cond[0]](self._s(cond[1])),self._s(cond[2])), p, evaluate=False)))
                        calculate_constant(cond, constant_det, constant_value, p, k, constant_lit)
                        if len(constant_det.keys()) == total_unknowns:
                            break                                           
                # finalize constants
                finalize_constants(constant_strings, reaction_strings, reactions, constant_value, constant_det)
                shear_block['constants'] = constant_det
                # plots
                plot_equations(plots, final_eqs, c, constant_value, constant_lit)
                shear_block['final_equations'] = final_eqs
                shear_block['plots'] = plots                          

            solution_blocks.append(shear_block)


        return solution_blocks

    def graph(self) -> go.Figure:
        """
        Creates a Plotly figure of the Beam.
        
        Returns:
            fig (go.Figure): The Plotly figure.
        """         
        BEAM_SIZE = self.beam_size
        ev = self._ev
        # scale down if too large
        scales = [10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
        for j, i in enumerate(scales):
            if self.beam_size > 2 * i and j == len(scales) - 1:
                BEAM_SIZE = self.beam_size/i
                ev = lambda l: self._ev(l)/i
                break       
            if self.beam_size > 2 * i and self.beam_size < 2 * scales[j + 1]:
                BEAM_SIZE = self.beam_size/i
                ev = lambda l: self._ev(l)/i
                break
        HEIGHT = (BEAM_HEIGHT * BEAM_SIZE)
        # Create the initial figure of the beam
        self.fig = plot_rectangle(BEAM_SIZE, HEIGHT, BEAM_COLOR)

        # Add the x and y axis arrows
        self.fig = add_axis_arrows(self.fig, HEIGHT)

        # Add points of interest
        total_points = len(list(self.points))
        for point in self.points:
            if point != 'A':
                test_pos = ev(self.points[chr(ord(point) - 1)])
            else:
                test_pos = 0            
            if test_pos < self.beam_size:
                self.fig = add_label(self.fig, ev(self.points[point]), -2 * HEIGHT, point, font_size = 16 if total_points < 6 else 12)
                if ord(point) > ord('A') and total_points > 2:
                    dist = str(simplify(self._s(f'({self.points[point]}) - ({self.points[chr(ord(point) - 1)]})')))
                    try:
                        dist = float(dist)
                        dist = f'{dist:1.2f}'
                    except:
                        if is_only_one_L_and_numbers(dist) and '/' not in dist:
                            dist = f'{ev(dist):1.2f}'
                    dist = dist.replace('.',',')
                    self.fig = add_hline_label(self.fig, -2 * HEIGHT, ev(self.points[chr(ord(point) - 1)]) + HEIGHT/8, ev(self.points[point]) - HEIGHT/8, format_label(dist).replace('*',''))

        # Add the links
        for link in self.links:
            if 'cantilever' in link:
                self.fig = add_cantilever(self.fig, HEIGHT, ev(link['cantilever']))
            if 'hinge' in link:
                self.fig = add_hinge(self.fig, HEIGHT, ev(link['hinge']))
            if 'mobile_support' in link:
                self.fig = add_mobile_support(self.fig, HEIGHT, ev(link['mobile_support']))
            if 'fixed_support' in link:
                self.fig = add_fixed_support(self.fig, HEIGHT, ev(link['fixed_support']))
            if 'roller' in link:
                self.fig = add_roller_support(self.fig, HEIGHT, ev(link['roller']))                
        
        # Add shear forces
        for force in self.shear_forces:
            if force['n'] < 0:
                if ev(force['start']) < self.beam_size:
                    if force['pos']:
                        if force['r']:
                            self.fig = add_vector(self.fig, (ev(force['start']), -1.5*HEIGHT + 0.15*HEIGHT), (ev(force['start']), -HEIGHT/2 + 0.15*HEIGHT), format_subs(force['value']), 'red')
                        else:
                            self.fig = add_vector(self.fig, (ev(force['start']), HEIGHT), (ev(force['start']), 2*HEIGHT), format_subs(force['value']))
                    else:
                        if force['r']:
                            self.fig = add_vector(self.fig, (ev(force['start']), HEIGHT/2), (ev(force['start']), -1.5*HEIGHT), format_subs(force['value']), 'red')
                        else:
                            self.fig = add_vector(self.fig, (ev(force['start']), 2.45*HEIGHT), (ev(force['start']), HEIGHT), format_subs(force['value']))
                else:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (ev(force['start']), 2.45*HEIGHT), (ev(force['start']), HEIGHT), format_subs(force['value']))                    
                    else:
                        self.fig = add_vector(self.fig, (ev(force['start']), HEIGHT), (ev(force['start']), 2*HEIGHT), format_subs(force['value']))
            else:
                x = np.linspace(ev(force['start']), ev(force['stop']), 1000)
                if force['n'] == 0:
                    y = [2* HEIGHT for i in range(len(x))]
                else:
                    f = lambdify(symbols('x'), self._get_load_equation(force).evalf(subs=self._get_symbol_value()), modules=[mapping, 'numpy'])
                    y = np.array(f(x))
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) * HEIGHT + HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, HEIGHT, format_subs(force['value']), up=force['pos'])
        
        # Add normal forces
        for force in self.normal_forces:
            if force['n'] < 0:
                if ev(force['start']) > 0:
                    if force['pos']:
                        if ev(force['start']) == BEAM_SIZE:
                            self.fig = add_vector(self.fig, (ev(force['start']), HEIGHT/2), (ev(force['start']) + HEIGHT, HEIGHT/2), format_subs(force['value']), 'orange')
                        else:
                            self.fig = add_vector(self.fig, (ev(force['start']) - HEIGHT, HEIGHT/2), (ev(force['start']), HEIGHT/2), format_subs(force['value']), 'orange')
                    else:
                        self.fig = add_vector(self.fig, (ev(force['start']) + HEIGHT, HEIGHT/2), (ev(force['start']), HEIGHT/2), format_subs(force['value']), 'orange')
                else:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (ev(force['start']), HEIGHT/2), (ev(force['start']) - HEIGHT, HEIGHT/2), format_subs(force['value']), 'orange')                 
                    else:
                        self.fig = add_vector(self.fig, (ev(force['start']) - HEIGHT, HEIGHT/2), (ev(force['start']), HEIGHT/2), format_subs(force['value']), 'orange')
            else:
                x = np.linspace(ev(force['start']), ev(force['stop']), 1000)
                if force['n'] == 0:
                    y = [2* HEIGHT for i in range(len(x))]
                else:
                    f = lambdify(symbols('x'), self._get_load_equation(force).evalf(subs=self._get_symbol_value()), modules=[mapping, 'numpy'])
                    y = np.array(f(x))
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) * HEIGHT + HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, HEIGHT, format_subs(force['value']), up=force['pos'], side=True)
        
        # Add twisting moments
        for moment in self.twisting_moments:
            if moment['n'] < 0:
                if ev(moment['start']) > 0:
                    if moment['pos']:
                        self.fig = add_vector(self.fig, (ev(moment['start']), HEIGHT/2), (ev(moment['start']) + 0.75*HEIGHT, HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (ev(moment['start']) + 0.7*HEIGHT, HEIGHT/2), (ev(moment['start']) + 1*HEIGHT, HEIGHT/2))
                    else:
                        self.fig = add_vector(self.fig, (ev(moment['start']) + 0.75*HEIGHT, HEIGHT/2), (ev(moment['start']), HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (ev(moment['start']) + 1*HEIGHT, HEIGHT/2), (ev(moment['start']) + 0.7*HEIGHT, HEIGHT/2))
                else:
                    if moment['pos']:
                        self.fig = add_vector(self.fig, (ev(moment['start']) + 0.75*HEIGHT, HEIGHT/2), (ev(moment['start']), HEIGHT/2), format_subs(moment['value'])) 
                        self.fig = add_vector(self.fig, (ev(moment['start']) + 1*HEIGHT, HEIGHT/2), (ev(moment['start']) + 0.7*HEIGHT, HEIGHT/2))                
                    else:
                        self.fig = add_vector(self.fig, (ev(moment['start']), HEIGHT/2), (ev(moment['start']) + 0.75*HEIGHT, HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (ev(moment['start']) + 0.7*HEIGHT, HEIGHT/2), (ev(moment['start']) + 1*HEIGHT, HEIGHT/2))
            else:
                x = np.linspace(ev(moment['start']), ev(moment['stop']), 1000)
                if moment['n'] == 0:
                    y = [2* HEIGHT for i in range(len(x))]
                else:
                    f = lambdify(symbols('x'), self._get_load_equation(moment).evalf(subs=self._get_symbol_value()), modules=[mapping, 'numpy'])
                    y = np.array(f(x))
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) * HEIGHT + HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, HEIGHT, format_subs(moment['value']), up=moment['pos'], side=True, double=True)    

        # Add bending moments
        for moment in self.bending_moments:
            if ev(moment['start']) < self.beam_size:
                if moment['pos']:
                    self.fig = add_semicircle_arrow(self.fig, ev(moment['start']), HEIGHT/2, orientation='clockwise', label=format_subs(moment['value']), radius=0.75*HEIGHT, color='rgba(255,140,0,1)')
                else:
                    self.fig = add_semicircle_arrow(self.fig, ev(moment['start']), HEIGHT/2, label=format_subs(moment['value']), radius=0.75*HEIGHT, color='rgba(255,140,0,1)')
            else:
                if moment['pos']:
                    self.fig = add_semicircle_arrow(self.fig, ev(moment['start']), HEIGHT/2, label=format_subs(moment['value']), radius=0.75*HEIGHT, color='rgba(255,140,0,1)')                 
                else:
                    self.fig = add_semicircle_arrow(self.fig, ev(moment['start']), HEIGHT/2, orientation='clockwise', label=format_subs(moment['value']), radius=0.75*HEIGHT, color='rgba(255,140,0,1)', start_angle=np.radians(90))

        return self.fig