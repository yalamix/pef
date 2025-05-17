from sqlalchemy import Column, Integer, String, ForeignKey, JSON, DateTime
from sqlalchemy.sql.expression import text
from sqlalchemy.orm import relationship
from .database import Base, engine
from datetime import datetime

class UserSession(Base):
    __tablename__ = "sessions"
    id = Column(String, primary_key=True, nullable=False)
    problem = relationship("Problem")

class Problem(Base):
    __tablename__ = "problems"
    id = Column(Integer, primary_key=True, nullable=False)
    parameters = Column(JSON, nullable=False)
    user_id = Column(String, ForeignKey("sessions.id", ondelete="CASCADE"))
    created_at = Column(DateTime,
                        nullable=False, server_default=text('CURRENT_TIMESTAMP'))

Base.metadata.create_all(engine)