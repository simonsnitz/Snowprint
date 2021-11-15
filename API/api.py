from api import app #, db
#from api.models import Sensor, Ligand, Reference


    #this is not working. Not sure why. Still need to manually import app, db when working in the command line
# @app.shell_context_processor
# def make_shell_context():
#     return {'db': db, 'Sensor': Sensor, 'Ligand': Ligand, 'Reference': Reference}

if __name__=="__main__":
    app.run(host="0.0.0.0", port=5000)