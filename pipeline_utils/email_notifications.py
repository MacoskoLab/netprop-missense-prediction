"""
Email notification utilities
"""

import smtplib
import ssl
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText


def send_email(body, subject, config):
    """Sends an email notification"""
    return  # Disabled for now

    # Load email settings from config
    try:
        sender_gmail = config["email"]["gmail_notification_email"]
        sender_gmail_password = config["email"]["gmail_notification_password"]
        receiver_email = config["email"]["email_notification_receiver"]
    except KeyError:
        print(
            "Some/all email settings not found in config. Skipping email notification."
        )
        return

    # Set up message
    message = MIMEMultipart()
    message["From"] = sender_gmail
    message["To"] = receiver_email
    message["Subject"] = subject
    message.attach(MIMEText(body, "plain"))

    # Set up SMTP server and send message
    context = ssl.create_default_context()
    with smtplib.SMTP_SSL("smtp.gmail.com", port=465, context=context) as server:
        server.login(sender_gmail, sender_gmail_password)
        server.sendmail(sender_gmail, receiver_email, message.as_string())
