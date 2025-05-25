from fastapi import FastAPI, Request, Response, Cookie, Depends, Form
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from sqlalchemy.orm import Session
from .database import get_db
from .interpreter import *
from .models import *
from .utils import *
from .graph import *
import logging
import uuid

app = FastAPI()

app.mount("/static", StaticFiles(directory="static"), name="static")

templates = Jinja2Templates(directory="templates")

@app.get("/", response_class=HTMLResponse)
async def index(response: Response, request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    if not session:
        session_id = str(uuid.uuid4())
        session = UserSession(
            id = session_id,
            problem = [Problem(
                parameters = DEFAULT_BEAM,
                user_id = session_id
            )]
        )
        db.add(session)
        db.commit()    
      
    BP = BeamProblem(json.dumps(session.problem[0].parameters))
    response = templates.TemplateResponse(
        request=request, name="index.html", response=response, context={'fig': fig_to_rotated_img(BP.graph())}, headers={'HX-Trigger':'renderCanvas'}
    )
    response.set_cookie(
        key="session_id",
        value=session_id,
        httponly=True,   # Not accessible via JavaScript.
        secure=False,    # Change to True in production (HTTPS).
        max_age=3600 * 720    # Session lifetime (e.g., 1 hour).
    )      
    return response

@app.get("/render", response_class=HTMLResponse)
async def render(request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()    
    BP = BeamProblem(json.dumps(session.problem[0].parameters))     
    return fig_to_rotated_img(BP.graph())

@app.get("/solve", response_class=HTMLResponse)
async def solve(request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()    
    BP = BeamProblem(json.dumps(session.problem[0].parameters))     
    result = BP.solve()
    return templates.TemplateResponse(
        request=request, name="solution.html", context={'result': result, 'len': len}
    )

@app.get("/get_variables", response_class=HTMLResponse)
async def get_variables(request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()  
    BP = BeamProblem(json.dumps(session.problem[0].parameters))     
    return templates.TemplateResponse(
        request=request, name="variables.html", context={'variables': BP.variables, 'format_subs': format_subs, 'latex': latex_with_threshold}
    )

@app.get("/get_links", response_class=HTMLResponse)
async def get_links(request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()  
    BP = BeamProblem(json.dumps(session.problem[0].parameters))     
    return templates.TemplateResponse(
        request=request, name="links.html", context={'links': BP.links, 'latex': latex_with_threshold}
    )

@app.get("/get_loads", response_class=HTMLResponse)
async def get_loads(request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()    
    BP = BeamProblem(json.dumps(session.problem[0].parameters))       
    return templates.TemplateResponse(
        request=request, name="loads.html", context={
            'shear_loads': BP.get_shear_forces(), 
            'normal_loads': BP.get_normal_forces(),
            'twisting_loads': BP.get_twisting_moments(),
            'bending_loads': BP.get_bending_moments()
        }
    )

@app.get("/get_conditions", response_class=HTMLResponse)
async def get_conditions(request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()    
    BP = BeamProblem(json.dumps(session.problem[0].parameters))       
    return templates.TemplateResponse(
        request=request, name="conditions.html", context={
            'conditions': BP.get_boundary_conditions()
        }
    )

@app.get("/debug", response_class=HTMLResponse)
async def debug(request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first() 
    return templates.TemplateResponse(
        request=request, name="debug.html", context={'msg': json.dumps(session.problem[0].parameters, indent=4)}
    )

@app.get("/reset")
async def reset(session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    session.problem[0].parameters = DEFAULT_BEAM
    db.commit()
    return RedirectResponse(url='/')

@app.get("/edit_variable/{variable}", response_class=HTMLResponse)
async def edit_variable(variable: str, request: Request, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    session = db.query(UserSession).where(UserSession.id == session_id).first()    
    BP = BeamProblem(json.dumps(session.problem[0].parameters))
    return templates.TemplateResponse(
        request=request, name="edit_variable.html", context={'variable':variable, 'value':BP.variables[variable]}
    )

# Modal to add elements

@app.get("/add_link", response_class=HTMLResponse)
async def add_link(request: Request):
    return templates.TemplateResponse(
        request=request, name="add_link.html",
    )

@app.get("/add_variable", response_class=HTMLResponse)
async def add_variable(request: Request):
    return templates.TemplateResponse(
        request=request, name="add_variable.html",
    )

@app.get("/add_shear", response_class=HTMLResponse)
async def add_shear(request: Request):
    return templates.TemplateResponse(
        request=request, name="add_force.html", context={'force': 'força cortante', 'route': '/shear'}
    )

@app.get("/add_normal", response_class=HTMLResponse)
async def add_normal(request: Request):
    return templates.TemplateResponse(
        request=request, name="add_force.html", context={'force': 'força normal', 'route': '/normal'}
    )

@app.get("/add_bending", response_class=HTMLResponse)
async def add_bending(request: Request):
    return templates.TemplateResponse(
        request=request, name="add_force.html", context={'force': 'momento fletor', 'route': '/bending'}
    )

@app.get("/add_twisting", response_class=HTMLResponse)
async def add_twisting(request: Request):
    return templates.TemplateResponse(
        request=request, name="add_force.html", context={'force': 'momento torsor', 'route': '/twisting'}
    )

# Post requests to add elements

@app.post("/link")
async def link(
        request: Request,
        link_type: str = Form(...),
        link_position: str = Form(...), 
        session_id: str = Cookie(None, alias="session_id"), 
        db: Session = Depends(get_db)
    ):
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    BP = BeamProblem(json.dumps(session.problem[0].parameters))
    result = BP.add_link(link_type, link_position)
    lt = {
        'cantilever': 'Engaste',
        'hinge': 'Rótula',
        'mobile_support': 'Apoio articulado fixo',
        'fixed_support': 'Apoio articulado móvel',
        'roller': 'Apoio móvel de roletes'
    }
    if result:
        session.problem[0].parameters = BP.to_dict()
        db.commit()
        context = {
            'msg': f'{lt[link_type]} adicionad{"o" if link_type != "hinge" else "a"}.',
            'success': result
        }
    else:
        context = {
            'msg': 'Erro ao adicionar vínculo. Confira se a posição está dentro do domínio ou se já não existe um vínculo nessa posição.',
            'success': result
        }        
    return templates.TemplateResponse(
        request=request, name="alert.html", context=context, headers={'HX-Trigger':'renderCanvas'}
    )

@app.post("/variable")
async def variable(
        request: Request,
        variable_name: str = Form(...),
        variable_value: float = Form(...),
        session_id: str = Cookie(None, alias="session_id"), 
        db: Session = Depends(get_db)        
    ):
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    BP = BeamProblem(json.dumps(session.problem[0].parameters))
    if variable_name not in BP.variables:
        context = {
            'msg': f'Adicionada variável {variable_name}. Valor: {variable_value}',
            'success': True
        } 
    else:
        context = {
            'msg': f'Variável {variable_name} atualizada. Valor: {variable_value}',
            'success': True
        }               
    BP.add_variable(variable_name, variable_value)
    session.problem[0].parameters = BP.to_dict()
    db.commit()
    return templates.TemplateResponse(
        request=request, name="alert.html", context=context, headers={'HX-Trigger':'renderCanvas'}
    )

@app.post("/shear", response_class=HTMLResponse)
async def shear(
        request: Request,
        force_value: str = Form(...),
        force_value_min: str = Form(...),
        pos: int = Form(...),
        n: int = Form(...),
        start: str = Form(...),
        stop: str = Form(...),
        session_id: str = Cookie(None, alias="session_id"), 
        db: Session = Depends(get_db)        
    ):
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    BP = BeamProblem(json.dumps(session.problem[0].parameters))
    BP.add_shear_force(force_value, force_value_min, start.replace(',','.'), stop.replace(',','.'), n, pos)
    session.problem[0].parameters = BP.to_dict()
    db.commit()    
    context = {
        'msg': f'Adicionada força cortante {"positiva" if pos else "negativa"}. Valor: {force_value}',
        'success': True
    }     
    return templates.TemplateResponse(
        request=request, name="alert.html", context=context, headers={'HX-Trigger':'renderCanvas'}
    )

@app.post("/normal", response_class=HTMLResponse)
async def normal(
        request: Request,
        force_value: str = Form(...),
        force_value_min: str = Form(...),
        pos: int = Form(...),
        n: int = Form(...),
        start: str = Form(...),
        stop: str = Form(...),
        session_id: str = Cookie(None, alias="session_id"), 
        db: Session = Depends(get_db)        
    ):
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    BP = BeamProblem(json.dumps(session.problem[0].parameters))
    BP.add_normal_force(force_value, force_value_min, start.replace(',','.'), stop.replace(',','.'), n, pos)
    session.problem[0].parameters = BP.to_dict()
    db.commit()    
    context = {
        'msg': f'Adicionada força normal {"positiva" if pos else "negativa"}. Valor: {force_value}',
        'success': True
    }     
    return templates.TemplateResponse(
        request=request, name="alert.html", context=context, headers={'HX-Trigger':'renderCanvas'}
    )

@app.post("/bending", response_class=HTMLResponse)
async def bending(
        request: Request,
        force_value: str = Form(...),
        force_value_min: str = Form(...),
        pos: int = Form(...),
        n: int = Form(...),
        start: str = Form(...),
        stop: str = Form(...),
        session_id: str = Cookie(None, alias="session_id"), 
        db: Session = Depends(get_db)        
    ):
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    BP = BeamProblem(json.dumps(session.problem[0].parameters))
    BP.add_bending_moment(force_value, force_value_min, start.replace(',','.'), stop.replace(',','.'), n, pos)
    session.problem[0].parameters = BP.to_dict()
    db.commit()    
    context = {
        'msg': f'Adicionada momento fletor {"positivo" if pos else "negativo"}. Valor: {force_value}',
        'success': True
    }     
    return templates.TemplateResponse(
        request=request, name="alert.html", context=context, headers={'HX-Trigger':'renderCanvas'}
    )

@app.post("/twisting", response_class=HTMLResponse)
async def twisting(
        request: Request,
        force_value: str = Form(...),
        force_value_min: str = Form(...),
        pos: int = Form(...),
        n: int = Form(...),
        start: str = Form(...),
        stop: str = Form(...),
        session_id: str = Cookie(None, alias="session_id"), 
        db: Session = Depends(get_db)        
    ):
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    BP = BeamProblem(json.dumps(session.problem[0].parameters))
    BP.add_twisting_moment(force_value, force_value_min, start.replace(',','.'), stop.replace(',','.'), n, pos)
    session.problem[0].parameters = BP.to_dict()
    db.commit()    
    context = {
        'msg': f'Adicionada momento torsor {"positivo" if pos else "negativo"}. Valor: {force_value}',
        'success': True
    }     
    return templates.TemplateResponse(
        request=request, name="alert.html", context=context, headers={'HX-Trigger':'renderCanvas'}
    )

@app.get("/read")
async def read_object(session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    # Retrieve the user's object from the session store.
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    if session:
        obj = json.dumps(session.problem.parameters)
    else:
        obj = {}  # Default object if session doesn't exist.
    return {"object": obj}

@app.post("/update")
async def update_object(response: Response, new_value: str, session_id: str = Cookie(None, alias="session_id"), db: Session = Depends(get_db)):
    # Create a new session if one doesn't exist.
    session = db.query(UserSession).where(UserSession.id == session_id).first()
    if not session:
        session_id = str(uuid.uuid4())
        session = UserSession(
            id = session_id,
            problem = [Problem(
                parameters = DEFAULT_BEAM,
                user_id = session_id
            )]
        )
        db.add(session)
        db.commit()
    
    # Update the object in the session store.
    session.problem = Problem(
                parameters = json.loads(new_value),
                user_id = session_id
            )

    # Set/update the session ID in the cookie.
    response.set_cookie(
        key="session_id",
        value=session_id,
        httponly=True,   # Not accessible via JavaScript.
        secure=False,    # Change to True in production (HTTPS).
        max_age=3600 * 720    # Session lifetime (e.g., 1 hour).
    )
    return {"message": "Object updated", "object": new_value}