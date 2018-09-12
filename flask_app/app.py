from flask import Flask
from flask import jsonify
from flask import request

app = Flask(__name__)

@app.route('/hello',methods=['POST'])
def hello():
    message = request.get_json(force=True)
    name = message['name']
    response = {
        'greeting': 'Hello, '+ name + '!'
    }
    
    return jsonify(response)
