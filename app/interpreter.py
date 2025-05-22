from typing import List
from .graph import *
from sympy import *
from .utils import *
import json

BEAM_HEIGHT = 0.5
BEAM_COLOR = 'blue'
DEFAULT_BEAM = {
    "beam_size": 5,
    "variables": {
        "L": 5
    },
    "links": [],
    "shear_forces": [],
    "normal_forces": [],
    "twisting_moments": [],
    "bending_moments": []
}

def format_subs(label: str) -> str:
    a = label.split('_')
    if len(a) > 1:
        total_subs = len(a) - 1
        label = '<sub>'.join(a)
        for i in range(total_subs):
            label += '</sub>'
    return label

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
        self.links = b['links']
        self.shear_forces = b['shear_forces']
        self.normal_forces = b['normal_forces']
        self.twisting_moments = b['twisting_moments']
        self.bending_moments = b['bending_moments']
    
    def to_dict(self) -> dict:
        """
        Returns a dictionary with the problem definitions.
        
        Returns:
            b: The problem dictionary.
        """              
        b = {}
        b['beam_size'] = self.beam_size
        b['variables'] = self.variables
        b['links'] = self.links
        b['shear_forces'] = self.shear_forces
        b['normal_forces'] = self.normal_forces
        b['twisting_moments'] = self.twisting_moments
        b['bending_moments'] = self.bending_moments
        return b
    
    def _update_variables(self, expr: str, default: float = 1) -> None:
        """
        Updates the variable list.   

        Parameters:
            expr        : Expression containing new variables.
            default        : Default value for new variables.
        """          
        new_symbols = parse_and_update_symbols(expr, list(self.variables.keys()))
        for s in new_symbols:
            if s not in self.variables:
                self.variables[s] = default     
    
    def _get_symbols(self) -> dict[Symbol]:
        """
        Converts the variables into Sympy Symbols.

        Returns:
            Symbol dictionary.
        """        
        return {name: Symbol(name) for name in list(self.variables.keys())}
        
    def _ev(self, expr: str) -> float:
        return float(sympify(expr, locals=self._get_symbols()).evalf(subs=self.variables))

    def assert_link_consistency(self, link_type: str, position: float) -> bool:
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
                self._update_variables(position)

        if self.assert_link_consistency(link_type, position):
            if link_type == 'cantilever':
                if self._ev(position) < self.beam_size/2:
                    position = 0
                else:
                    position = self.beam_size
            if self._ev(position) < 0:
                position = 0
            if self._ev(position) > self.beam_size:
                position = self.beam_size
            self.links.append({link_type: position})
            return True
        return False

    def _create_load(self, value: str, start: str, stop: str, n: int, pos: bool = True) -> dict:
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
                self._update_variables(value)
        if start not in self.variables:
            try:
                float(start)
            except:     
                self._update_variables(start)
        if stop not in self.variables:
            try:
                float(stop)
            except:     
                self._update_variables(stop, 2)

        # if n >= 0 assert stop > start
        if n >= 0:
            if self._ev(stop) < self._ev(start):
                self.stop = 'L'
        # constrain load to problem domain
        if self._ev(start) < 0:
            start = 0
        if self._ev(stop) > self.beam_size:
            stop = 'L'

        return {
                "value": value,
                "start": start,
                "stop": stop,
                "n": n,
                "pos": pos
            }          

    def add_shear_force(self, value: str, start: str, stop: str, n: int, pos: bool = True):
        """
        Adds a concentrated or distributed shear force.   

        Parameters:
            value        : Value of the force.        
            start         : Start position of the domain of the force.
            stop         : End position of the domain of the force.
            n         : Exponent of the polynomial representing the force.
            pos         : True if the force is positive, False otherwise.
        """ 
        self.shear_forces.append(self._create_load(value, start, stop, n, pos))  

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
    
    def _get_load_expression(self, force: dict) -> str:
        """
        Generates the latex string of the singularity function representing the load.   

        Parameters:
            force        : Load dictionary.        

        Returns:
            out: Latex string of the load expression.
        """          
        x = symbols('x', positive=True)
        try:
            p = float(force['value'])
        except:
            p = sympify(force['value'], locals=self._get_symbols())            
        if self._ev(force['stop']) < self.beam_size and force['n'] >= 0:
            eq = p * SingularityFunction(x, force['start'], force['n']) - p * SingularityFunction(x, force['stop'], force['n'])
        else:
            eq = p * SingularityFunction(x, force['start'], force['n'])
        return latex_with_threshold(eq)
    
    def _load_expressions(self, loads: List[dict]) -> List[str]:
        """
        Returns a list of the latex form of the loads.
        
        Returns:
            b: The list of latex strings.
        """        
        return [self._get_load_expression(load) for load in loads]        

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
        ...

    def graph(self) -> go.Figure:
        """
        Creates a Plotly figure of the Beam.
        
        Returns:
            fig (go.Figure): The Plotly figure.
        """         
        # Create the initial figure of the beam
        self.fig = plot_rectangle(self.beam_size, BEAM_HEIGHT, BEAM_COLOR)

        # Add the x and y axis arrows
        self.fig = add_axis_arrows(self.fig)

        # Add the links
        for link in self.links:
            if 'cantilever' in link:
                self.fig = add_cantilever(self.fig, BEAM_HEIGHT, self._ev(link['cantilever']))
            if 'hinge' in link:
                self.fig = add_hinge(self.fig, BEAM_HEIGHT, self._ev(link['hinge']))
            if 'mobile_support' in link:
                self.fig = add_mobile_support(self.fig, BEAM_HEIGHT, self._ev(link['mobile_support']))
            if 'fixed_support' in link:
                self.fig = add_fixed_support(self.fig, BEAM_HEIGHT, self._ev(link['fixed_support']))
            if 'roller' in link:
                self.fig = add_roller_support(self.fig, BEAM_HEIGHT, self._ev(link['roller']))                
        
        # Add shear forces
        for force in self.shear_forces:
            if force['n'] < 0:
                if self._ev(force['start']) < self.beam_size:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), BEAM_HEIGHT), (self._ev(force['start']), BEAM_HEIGHT + 1), format_subs(force['value']))
                    else:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), BEAM_HEIGHT + 1), (self._ev(force['start']), BEAM_HEIGHT), format_subs(force['value']))
                else:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), BEAM_HEIGHT + 1), (self._ev(force['start']), BEAM_HEIGHT), format_subs(force['value']))                    
                    else:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), BEAM_HEIGHT), (self._ev(force['start']), BEAM_HEIGHT + 1), format_subs(force['value']))
            else:
                x = np.linspace(self._ev(force['start']), self._ev(force['stop']), 100)
                if force['n'] == 0:
                    y = [2* BEAM_HEIGHT for i in range(len(x))]
                else:
                    y = (x - self._ev(force['start'])) ** force['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + BEAM_HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, BEAM_HEIGHT, format_subs(force['value']), up=force['pos'])
        
        # Add normal forces
        for force in self.normal_forces:
            if force['n'] < 0:
                if self._ev(force['start']) > 0:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), BEAM_HEIGHT/2), (self._ev(force['start']) + 1, BEAM_HEIGHT/2), format_subs(force['value']))
                    else:
                        self.fig = add_vector(self.fig, (self._ev(force['start']) + 1, BEAM_HEIGHT/2), (self._ev(force['start']), BEAM_HEIGHT/2), format_subs(force['value']))
                else:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (self._ev(force['start']) + 1, BEAM_HEIGHT/2), (self._ev(force['start']), BEAM_HEIGHT/2), format_subs(force['value']))                 
                    else:
                        self.fig = add_vector(self.fig, (self._ev(force['start']), BEAM_HEIGHT/2), (self._ev(force['start']) + 1, BEAM_HEIGHT/2), format_subs(force['value']))
            else:
                x = np.linspace(self._ev(force['start']), self._ev(force['stop']), 100)
                if force['n'] == 0:
                    y = [2* BEAM_HEIGHT for i in range(len(x))]
                else:
                    y = (x - self._ev(force['start'])) ** force['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + BEAM_HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, BEAM_HEIGHT, format_subs(force['value']), up=force['pos'], side=True)
        
        # Add twisting moments
        for moment in self.twisting_moments:
            if moment['n'] < 0:
                if self._ev(moment['start']) > 0:
                    if moment['pos']:
                        self.fig = add_vector(self.fig, (self._ev(moment['start']), BEAM_HEIGHT/2), (self._ev(moment['start']) + 1, BEAM_HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 0.5, BEAM_HEIGHT/2), (self._ev(moment['start']) + 1.5, BEAM_HEIGHT/2), format_subs(moment['value']))
                    else:
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 1, BEAM_HEIGHT/2), (self._ev(moment['start']), BEAM_HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 1.5, BEAM_HEIGHT/2), (self._ev(moment['start']) + 0.5, BEAM_HEIGHT/2), format_subs(moment['value']))
                else:
                    if moment['pos']:
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 1, BEAM_HEIGHT/2), (self._ev(moment['start']), BEAM_HEIGHT/2), format_subs(moment['value'])) 
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 1.5, BEAM_HEIGHT/2), (self._ev(moment['start']) + 0.5, BEAM_HEIGHT/2), format_subs(moment['value']))                
                    else:
                        self.fig = add_vector(self.fig, (self._ev(moment['start']), BEAM_HEIGHT/2), (self._ev(moment['start']) + 1, BEAM_HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (self._ev(moment['start']) + 0.5, BEAM_HEIGHT/2), (self._ev(moment['start']) + 1.5, BEAM_HEIGHT/2), format_subs(moment['value']))
            else:
                x = np.linspace(self._ev(moment['start']), self._ev(moment['stop']), 100)
                if moment['n'] == 0:
                    y = [2* BEAM_HEIGHT for i in range(len(x))]
                else:
                    y = (x - self._ev(moment['start'])) ** moment['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + BEAM_HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, BEAM_HEIGHT, format_subs(moment['value']), up=moment['pos'], side=True, double=True)    

        # Add bending moments
        for moment in self.bending_moments:
            if self._ev(moment['start']) < self.beam_size:
                if moment['pos']:
                    self.fig = add_semicircle_arrow(self.fig, self._ev(moment['start']), BEAM_HEIGHT/2, orientation='clockwise', label=format_subs(moment['value']))
                else:
                    self.fig = add_semicircle_arrow(self.fig, self._ev(moment['start']), BEAM_HEIGHT/2, label=format_subs(moment['value']))
            else:
                if moment['pos']:
                    self.fig = add_semicircle_arrow(self.fig, self._ev(moment['start']), BEAM_HEIGHT/2, label=format_subs(moment['value']))                 
                else:
                    self.fig = add_semicircle_arrow(self.fig, self._ev(moment['start']), BEAM_HEIGHT/2, orientation='clockwise', label=format_subs(moment['value']))

        return self.fig