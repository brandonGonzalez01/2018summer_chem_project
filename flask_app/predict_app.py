
from flask import request, jsonify, Flask

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit import DataStructs
from rdkit.Chem import AllChem

import numpy as np
import io
import rdkit
import keras

from keras.models import Sequential
from keras.models import load_model

import pubchempy as pcp

app = Flask(__name__)
app.debug = False

def get_model():
    global model
    global excipientModel
    global model_PH
    excipientModel = load_model('chem_E_model.h5')
    model = load_model('chem_logs_model.h5')
    model_PH=load_model('chem_E_model.h5')
    x = preprocess_SMILE("C")
    model.predict(x)
    model_PH.predict(x)
    excipientModel.predict(x)    
    print(" * Model loaded!")

def get_weight(smiles):
    #check if smiles was single SMILES input
    if(type(smiles) is str):
        return Chem.Descriptors.MolWt(Chem.MolFromSmiles(smiles))
    #check if SMILES smiles was from Excipient option
    elif(smiles[0]=='E'):
        return Chem.Descriptors.MolWt(Chem.MolFromSmiles(smiles[1]))
    else:
        #else smiles is an array of smiles
        weights = []
        for i in range(len(smiles)):
            mol = Chem.MolFromSmiles(smiles[i])
            wieght = Chem.Descriptors.MolWt(mol)
            weights.append(wieght)
        return weights

def get_name(smiles):
    if(type(smiles) is str):
        try:
            syn = pcp.get_synonyms(smiles,'smiles')[0]['Synonym']
            print(syn)
            return syn[0]
        except:
            try:
                return pcp.get_properties('IUPACName',smiles,'smiles')[0]['IUPACName']
            except:
                return smiles
    else:
        names = []
        for i in range(len(smiles)):
            if(smiles[i]==''):
                continue
            else:
                try:
                    syn = pcp.get_synonyms(smiles[i],'smiles')[0]['Synonym'][0]
                    print('syn @ ',i,' : ',syn)
                    names.append(syn)
                except:
                    try:
                        IUPAC = pcp.get_properties('IUPACName',smiles[i],'smiles')
                        print(f"{i}-{smiles[i]} , IUPAC used instead {IUPAC}" )
                        names.append(IUPAC[0]['IUPACName'][0])
                    except:
                        names.append(smiles[i][:10])
    return names

def preprocess_SMILE(smiles):
    print(smiles)
    fpArray = []
    #length = number of bits the Fingerprint has
    length = 2048
    #print('length: ',len(smiles),smiles)
    if(type(smiles) is str or smiles[0]=="E" ):
        mol = ""
        
        #check if excipient option is selected
        if smiles[0]=="E":
            mol = Chem.MolFromSmiles(smiles[1])
        else:
            mol = Chem.MolFromSmiles(smiles)
        #error check
        if(mol is None):
            print('mol invalid: ',smiles)
        else:
            fp = AllChem.rdmolops.RDKFingerprint(mol,fpSize = length)
            fpArray.append(fp)  
    
    else:#else the smiles in an array for the water solubility model
        #this split function bc the first array SMILE has an abberanr " in it 
        print("length of smiles: ",len(smiles))
        for i in range(len(smiles)):
            mol = Chem.MolFromSmiles(smiles[i])
            if(mol is None):
                continue
            else:
                fp = AllChem.rdmolops.RDKFingerprint(mol, fpSize = length)
                fpArray.append(fp)
    #convert the fingerprint array into a numpy array
    X = np.array(fpArray) 
    print(X)
    return X

print("* Loading model...")


get_model()

@app.route("/predict",methods=['POST'])

def predict():
    response={}
    message = request.get_json(force=True)
    print(message)
    #check if the chosen option was to use the Excipient model
    keys = list(message.keys())
    if('toBeNamed' in keys ):
        print(message['toBeNamed'])
        print('success')
        encoded_smile=message['toBeNamed']
        name = get_name(encoded_smile)
        response = {
            'names': name
        }
        return jsonify(response)
    #check if the chosen option was to use the Excipient model message['mol'][0]=="E" is true
    elif(message['mol'][0] and message['mol'][0]=="E"):
        print(message['mol'][0])
        encoded_smile= message['mol'][1]
        
        #get weight properties of the SMILES
        weight = get_weight(encoded_smile)
        
        #get the name of the SMILES in iupac
        name = get_name(encoded_smile)
        
        processed_smile = preprocess_SMILE(encoded_smile)
        prediction = excipientModel.predict(processed_smile)
        print(prediction.tolist())
        response={
            'prediction': prediction.tolist(),
            'weight' : weight,
            'name' : name
        }
        return jsonify(response)
    else:     
        
        #else this model was the solubility in water model
        encoded_smile = message['mol']
        
        #get the weights
        weight = get_weight(encoded_smile)
        
        #get the names of the encoded_smile smiles list
        name = get_name(encoded_smile)

        processed_smile = preprocess_SMILE(encoded_smile)
        prediction = model.predict(processed_smile)
        print(prediction.tolist())
        response = {
                'prediction' : prediction.tolist(),     
                'weight' : weight,
                'name' : name,
            }
        print(response)
        return jsonify(response)

@app.route("/pH",methods=['POST'])
#print("pH page accesed")
def ph_page():
    return 0


        
    

    