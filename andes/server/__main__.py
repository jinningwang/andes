import os
import andes
import threading  # NOQA

from andes.utils.math import to_number
import flask
from flask import Flask, request  # NOQA
from flask_restful import Resource, Api  # NOQA
from flask.json import jsonify  # NOQA

app = Flask(__name__)
api = Api(app)


andes.main.config_logger()
systems = {}
sim_thread = {}


@app.route('/status')
def get_status():
    sysid = request.args.get('sysid', None)

    if sysid not in systems:
        return jsonify('0')
    else:
        return jsonify('1')


@app.route('/load')
def load():
    default_path = os.getcwd()
    name = request.args.get('name', '')
    path = os.path.join(default_path, name)

    n_system = len(systems)
    try:
        system_instance = andes.main.run(case=path)
        globals()['systems'][str(n_system + 1)] = system_instance

    except FileNotFoundError:
        flask.abort(404)

    return jsonify(len(systems))


@app.route('/run')
def run():
    sysid = request.args.get('sysid', None)
    simulation_time = request.args.get('time', 0)

    if sysid not in systems:
        flask.abort(400)
    if simulation_time == 0:
        flask.abort(400)

    system = systems[sysid]

    if request.method == "GET":
        system.pflow.run()
        system.tds.init()

        system.tds.config.qrt = True
        system.tds.config.tf = float(simulation_time)

        thread = threading.Thread(target=system.tds.run)
        sim_thread[sysid] = thread
        thread.start()

    return jsonify({'response': 'success'})


@app.route('/param', methods=['GET', 'POST'])
def get_model_param():

    sysid = request.args.get('sysid', None)
    model_name = request.args.get('name', None)
    var_name = request.args.get('var', None)
    idx = request.args.get('idx', None)
    value = request.args.get('value', None)
    sysbase = request.args.get('sysbase', 'False')

    if not sysid:
        flask.abort(400)
    if sysbase == 'False':
        sysbase = False
    elif sysbase == 'True':
        sysbase = True

    system = systems[sysid]

    if request.method == 'GET':

        if idx is not None:
            idx = to_number(idx)

        if request.method == 'GET':
            if not model_name or (model_name not in system.devman.devices):
                return 'Model name <{}> invalid or not loaded in system'.format(model_name)

            else:  # with model_name

                model_ref = system.__dict__[model_name]
                if var_name:  # with `var_name`
                    if var_name not in model_ref.__dict__:
                        return 'Error: variable <{}> not exist in <{}>'.format(var_name, model_name)

                    if idx:
                        return jsonify(model_ref.get_field(field=var_name, idx=idx))
                    elif not idx:
                        return jsonify(list(model_ref.get_field(field=var_name)))

                elif not var_name:  # without `var_name`
                    if idx:
                        if idx not in model_ref.idx:
                            return 'Error: idx <{}> not exist in <{}>'.format(idx, model_name)
                        return jsonify(model_ref.get_element_data(idx))
                    elif not idx:
                        return jsonify(model_ref.data_to_list())

    elif request.method == 'POST':
        if any([model_name, var_name, idx, value]) is None:
            flask.abort(400)

        model_ref = system.__dict__[model_name]

        if var_name not in model_ref.__dict__:
            flask.abort(404)

        idx = to_number(idx)
        if idx not in model_ref.idx:
            flask.abort(404)

        if sysbase == 'False':
            sysbase = False
        elif sysbase == 'True':
            sysbase = True

        model_ref.set_field(var_name, idx, value, sysbase)
        model_ref.reload_new_param()
        return jsonify(model_ref.get_field(var_name, idx))


@app.route('/sim_time')
def get_simulation_time():
    """
    Get the simulation time of a system

    Returns
    -------
    str : time
    """
    sysid = request.args.get('sysid', None)
    if sysid is None:
        flask.abort(400)

    system = systems['sysid']

    return jsonify(system.tds.t)


@app.route('/streaming')
def get_streaming_data():
    """
    Get the latest variable in `varout`

    Returns
    -------

    """
    sysid = request.args.get('sysid', None)
    if sysid is None:
        flask.abort(400)

    system = systems[sysid]

    return jsonify(system.varout.get_latest_data())


if __name__ == '__main__':
    app.run(port='5000')
