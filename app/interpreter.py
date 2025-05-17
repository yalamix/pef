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
          
    def add_link(self, link_type: str, position: float) -> bool:
        """
        Adds a new link.

        Parameters:
            link_type        : Name of the link.        
            position         : Position of the link.             
        
        Returns:
            True if the link was added, False otherwise.
        """        
        if self.assert_link_consistency(link_type, position):
            self.links.append({link_type: position})
            return True
        return False
    
    def add_shear_force(self, value: str, start: float, stop: float, n: int, pos: bool = True):
        """
        Adds a concentrated or distributed shear force.   

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
                self.variables[value] = 1
        self.shear_forces.append(
            {
                "value": value,
                "start": start,
                "stop": stop,
                "n": n,
                "pos": pos
            }            
        )  

    def add_normal_force(self, value: str, start: float, stop: float, n: int, pos: bool = True):
        """
        Adds a concentrated or distributed normal force.   

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
                self.variables[value] = 1        
        self.normal_forces.append(
            {
                "value": value,
                "start": start,
                "stop": stop,
                "n": n,
                "pos": pos
            }            
        )  

    def add_twisting_moment(self, value: str, start: float, stop: float, n: int, pos: bool = True):
        """
        Adds a concentrated or distributed twisting moment.   

        Parameters:
            value        : Value of the moment.        
            start         : Start position of the domain of the moment.
            stop         : End position of the domain of the moment.
            n         : Exponent of the polynomial representing the moment.
            pos         : True if the moment is positive, False otherwise.
        """ 
        if value not in self.variables:
            try:
                float(value)
            except:
                self.variables[value] = 1        
        self.twisting_moments.append(
            {
                "value": value,
                "start": start,
                "stop": stop,
                "n": n,
                "pos": pos
            }            
        )  

    def add_bending_moment(self, value: str, start: float, stop: float, n: int, pos: bool = True):
        """
        Adds a concentrated or distributed bending moment.   

        Parameters:
            value        : Value of the moment.        
            start         : Start position of the domain of the moment.
            stop         : End position of the domain of the moment.
            n         : Exponent of the polynomial representing the moment.
            pos         : True if the moment is positive, False otherwise.
        """ 
        if value not in self.variables:
            try:
                float(value)
            except:
                self.variables[value] = 1        
        self.bending_moments.append(
            {
                "value": value,
                "start": start,
                "stop": stop,
                "n": n,
                "pos": pos
            }            
        )  
    
    def get_normal_forces(self) -> List[str]:
        """
        Returns a list of the latex form of the normal forces.
        
        Returns:
            b: The list of latex strings.
        """
        forces = []
        x = symbols('x', positive=True)
        for force in self.normal_forces:
            try:
                p = float(force['value'])
            except:
                p = symbols(force['value'])            
            if force['stop'] < self.beam_size and force['n'] >= 0:
                eq = p * SingularityFunction(x, force['start'], force['n']) - p * SingularityFunction(x, force['stop'], force['n'])
            else:
                eq = p * SingularityFunction(x, force['start'], force['n'])
            forces.append(latex_with_threshold(eq))
        return forces

    def get_shear_forces(self) -> List[str]:
        """
        Returns a list of the latex form of the shear forces.
        
        Returns:
            b: The list of latex strings.
        """
        forces = []
        x = symbols('x', positive=True)
        for force in self.shear_forces:
            try:
                p = float(force['value'])
            except:
                p = symbols(force['value'])            
            if force['stop'] < self.beam_size and force['n'] >= 0:
                eq = p * SingularityFunction(x, force['start'], force['n']) - p * SingularityFunction(x, force['stop'], force['n'])
            else:
                eq = p * SingularityFunction(x, force['start'], force['n'])
            forces.append(latex_with_threshold(eq))
        return forces    

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
                self.fig = add_cantilever(self.fig, BEAM_HEIGHT, link['cantilever'])
            if 'hinge' in link:
                self.fig = add_hinge(self.fig, BEAM_HEIGHT, link['hinge'])
            if 'mobile_support' in link:
                self.fig = add_mobile_support(self.fig, BEAM_HEIGHT, link['mobile_support'])
            if 'fixed_support' in link:
                self.fig = add_fixed_support(self.fig, BEAM_HEIGHT, link['fixed_support'])
            if 'roller' in link:
                self.fig = add_roller_support(self.fig, BEAM_HEIGHT, link['roller'])                
        
        # Add shear forces
        for force in self.shear_forces:
            if force['n'] < 0:
                if force['start'] < self.beam_size:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (force['start'], BEAM_HEIGHT), (force['start'], BEAM_HEIGHT + 1), format_subs(force['value']))
                    else:
                        self.fig = add_vector(self.fig, (force['start'], BEAM_HEIGHT + 1), (force['start'], BEAM_HEIGHT), format_subs(force['value']))
                else:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (force['start'], BEAM_HEIGHT + 1), (force['start'], BEAM_HEIGHT), format_subs(force['value']))                    
                    else:
                        self.fig = add_vector(self.fig, (force['start'], BEAM_HEIGHT), (force['start'], BEAM_HEIGHT + 1), format_subs(force['value']))
            else:
                x = np.linspace(force['start'], force['stop'], 100)
                if force['n'] == 0:
                    y = [2* BEAM_HEIGHT for i in range(len(x))]
                else:
                    y = (x - force['start']) ** force['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + BEAM_HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, BEAM_HEIGHT, format_subs(force['value']), up=force['pos'])
        
        # Add normal forces
        for force in self.normal_forces:
            if force['n'] < 0:
                if force['start'] > 0:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (force['start'], BEAM_HEIGHT/2), (force['start'] + 1, BEAM_HEIGHT/2), format_subs(force['value']))
                    else:
                        self.fig = add_vector(self.fig, (force['start'] + 1, BEAM_HEIGHT/2), (force['start'], BEAM_HEIGHT/2), format_subs(force['value']))
                else:
                    if force['pos']:
                        self.fig = add_vector(self.fig, (force['start'] + 1, BEAM_HEIGHT/2), (force['start'], BEAM_HEIGHT/2), format_subs(force['value']))                 
                    else:
                        self.fig = add_vector(self.fig, (force['start'], BEAM_HEIGHT/2), (force['start'] + 1, BEAM_HEIGHT/2), format_subs(force['value']))
            else:
                x = np.linspace(force['start'], force['stop'], 100)
                if force['n'] == 0:
                    y = [2* BEAM_HEIGHT for i in range(len(x))]
                else:
                    y = (x - force['start']) ** force['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + BEAM_HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, BEAM_HEIGHT, format_subs(force['value']), up=force['pos'], side=True)
        
        # Add twisting moments
        for moment in self.twisting_moments:
            if moment['n'] < 0:
                if moment['start'] > 0:
                    if moment['pos']:
                        self.fig = add_vector(self.fig, (moment['start'], BEAM_HEIGHT/2), (moment['start'] + 1, BEAM_HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (moment['start'] + 0.5, BEAM_HEIGHT/2), (moment['start'] + 1.5, BEAM_HEIGHT/2), format_subs(moment['value']))
                    else:
                        self.fig = add_vector(self.fig, (moment['start'] + 1, BEAM_HEIGHT/2), (moment['start'], BEAM_HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (moment['start'] + 1.5, BEAM_HEIGHT/2), (moment['start'] + 0.5, BEAM_HEIGHT/2), format_subs(moment['value']))
                else:
                    if moment['pos']:
                        self.fig = add_vector(self.fig, (moment['start'] + 1, BEAM_HEIGHT/2), (moment['start'], BEAM_HEIGHT/2), format_subs(moment['value'])) 
                        self.fig = add_vector(self.fig, (moment['start'] + 1.5, BEAM_HEIGHT/2), (moment['start'] + 0.5, BEAM_HEIGHT/2), format_subs(moment['value']))                
                    else:
                        self.fig = add_vector(self.fig, (moment['start'], BEAM_HEIGHT/2), (moment['start'] + 1, BEAM_HEIGHT/2), format_subs(moment['value']))
                        self.fig = add_vector(self.fig, (moment['start'] + 0.5, BEAM_HEIGHT/2), (moment['start'] + 1.5, BEAM_HEIGHT/2), format_subs(moment['value']))
            else:
                x = np.linspace(moment['start'], moment['stop'], 100)
                if moment['n'] == 0:
                    y = [2* BEAM_HEIGHT for i in range(len(x))]
                else:
                    y = (x - moment['start']) ** moment['n']
                    y = (y-np.min(y))/(np.max(y)-np.min(y)) + BEAM_HEIGHT
                self.fig = add_curve_with_y_cutoff_fill(self.fig, x, y, BEAM_HEIGHT, format_subs(moment['value']), up=moment['pos'], side=True, double=True)    

        # Add bending moments
        for moment in self.bending_moments:
            if moment['start'] < self.beam_size:
                if moment['pos']:
                    self.fig = add_semicircle_arrow(self.fig, moment['start'], BEAM_HEIGHT/2, orientation='clockwise', label=format_subs(moment['value']))
                else:
                    self.fig = add_semicircle_arrow(self.fig, moment['start'], BEAM_HEIGHT/2, label=format_subs(moment['value']))
            else:
                if moment['pos']:
                    self.fig = add_semicircle_arrow(self.fig, moment['start'], BEAM_HEIGHT/2, label=format_subs(moment['value']))                 
                else:
                    self.fig = add_semicircle_arrow(self.fig, moment['start'], BEAM_HEIGHT/2, orientation='clockwise', label=format_subs(moment['value']))

        return self.fig