from ASB_app import app, scheduler

if __name__ == '__main__':
    scheduler.start()
    app.run(debug=True, host='0.0.0.0', port=5000, use_reloader=False)
