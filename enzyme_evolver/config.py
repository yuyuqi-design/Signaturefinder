import os
basedir = os.path.abspath(os.path.dirname(__file__))
import redis

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'testing_key'
    SECURITY_PASSWORD_SALT = os.environ.get("SECURITY_PASSWORD_SALT", '8439842')

    PRODUCTION = os.environ.get('PRODUCTION') or False

    REDIS_URL = os.environ.get('REDIS_URL') or 'redis://'
    MONGODB_HOST = os.environ.get('MONGO_HOST') or 'localhost'
    MONGODB_DB = os.environ.get('MONGODB_DB') or 'enzyme_evolver_db'
    MONGODB_PORT = os.environ.get('MONGODB_PORT') or 27017

    ADMIN_EMAIL = os.environ.get('ADMIN_EMAIL') or 'admin@email.com'
    ADMIN_PASSWORD = os.environ.get('ADMIN_PASSWORD') or 'password'

    USE_EMAIL_CONFIRMATION = os.environ.get('USE_EMAIL_CONFIRMATION') or False
    SECURITY_REGISTERABLE = True
    SECURITY_SEND_REGISTER_EMAIL = USE_EMAIL_CONFIRMATION
    SECURITY_CONFIRMABLE = USE_EMAIL_CONFIRMATION
    SECURITY_CHANGEABLE = USE_EMAIL_CONFIRMATION or True
    SECURITY_RECOVERABLE = USE_EMAIL_CONFIRMATION
    SECURITY_LOGIN_WITHOUT_CONFIRMATION = not USE_EMAIL_CONFIRMATION
    SECURITY_CONFIRMABLE = True
    SECURITY_RECOVERABLE =True

    if SECURITY_CONFIRMABLE == True:
        SECURITY_POST_REGISTER_VIEW = '/confirm'
        SECURITY_POST_CHANGE_VIEW = '/change'
    else:
        SECURITY_POST_REGISTER_VIEW = '/'
        SECURITY_POST_CHANGE_VIEW = '/'


    MAIL_SERVER = os.environ.get('MAIL_SERVER') or 'smtp-mail.outlook.com'
    MAIL_PORT = os.environ.get('MAIL_PORT') or '587'
    MAIL_USE_SSL = False
    MAIL_USE_TLS = True
    MAIL_DEBUG = True
    MAIL_USERNAME = os.environ.get('EMAIL_ADDRESS') or 'enzyme_evolver@outlook.com'
    MAIL_PASSWORD = os.environ.get('EMAIL_PASSWORD') or 'enzymeevolverMIB'
    MAIL_DEFAULT_SENDER = 'enzyme_evolver@outlook.com'
    MAIL_MAX_EMAILS = None
    #MAIL_SUPPRESS_SEND = False
    MAIL_ASCII_ATTACHMENTS = False

    CACHE_TYPE = "redis"
    CACHE_DEFAULT_TIMEOUT = 300
    CACHE_REDIS_HOST = REDIS_URL

    SESSION_TYPE = "redis"
    SESSION_REDIS = redis.from_url(REDIS_URL)
    SESSION_USE_SIGNER = True

